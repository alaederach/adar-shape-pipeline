#!/usr/bin/env python3

"""
Given a modified-condition and untreated-condition coordinate-sorted BAM,
identify QNAMEs whose molecules contain the strongest possible
refbase->readbase mutation pattern while retaining at least a minimum read
depth.

Supports multiple genomic regions in a single run via numbered argument pairs
(--start1/--end1/--min-rate1, --start2/--end2/--min-rate2, etc.). Each region
is processed independently.

Candidate positions for each region may be derived automatically from the
modified-condition profile.txt using a per-region --min-rateN threshold, or
supplied directly as a comma-separated list of 1-indexed positions via
--candidatesN. --min-rateN and --candidatesN are mutually exclusive per region.

For each region:
  - Candidate positions are selected either from the modified-condition profile.txt
    (profile + min-rate mode) or from user-supplied positions (--candidatesN mode).
  - A --minimum-coverage pre-filter restricts analysis to reads spanning at least
    the specified percentage of the region (applied to all counting steps).
  - The threshold is determined jointly from both BAMs (highest value where both
    exceed --edited-min-depth).
  - Two modified-condition BAM output files are written (one per input BAM).
  - Optionally (default on, disable with --no-unedited), two unmodified BAM
    output files are written. A molecule qualifies as unmodified if no mate has
    a mutation at a candidate position and the union of non-candidate mutations
    across all mates does not exceed --unmodified-tolerance.

A single log file capturing run parameters and per-region summaries is written
to --out-dir, named <mod_stem>_log.txt where <mod_stem> is the modbam filename
with all dot-suffixes stripped (e.g. sample.unique.sorted.bam -> sample).

Read counting and mutation scoring
-----------------------------------
Molecules are scored at the QNAME level, not the read level. BAMs may contain
a mixture of merged single-end reads and unmerged concordant pairs sharing the
same QNAME. For any given QNAME, all mates that pass the minimum-coverage
filter contribute to the molecule's score:

  - Edited scoring: count of distinct candidate positions showing readbase
    across the union of all mates.
  - Unmodified scoring: count of distinct refbase positions covered across the
    union of all mates. A candidate-position mutation on any mate disqualifies
    the entire molecule. Non-candidate mutations are also pooled across mates
    and evaluated against --unmodified-tolerance as a union.

Coverage qualification is evaluated at the QNAME level using the
non-overlapping union of all mates' reference spans within the region. A QNAME
qualifies if that union meets or exceeds --minimum-coverage % of the region
length. This correctly handles the case where neither mate individually spans
enough of the region but together they do (e.g., R1 covers the left portion,
R2 the right, with a gap in between). All mates of a qualified QNAME
contribute positions regardless of their individual span. QNAMEs where no mate
overlaps the region, or whose union span falls below the threshold, are
excluded from the output.

Assumptions
-----------
1. Both BAMs contain only uniquely mapped, primary-mapped reads (no secondary
   or supplementary alignments). Reads may be single-end (merged) or unmerged
   concordant pairs; both are handled correctly.
2. If --modprofile is provided, it must be tab-delimited with columns: Nucleotide,
   Sequence, Modified_rate. Required unless all regions use --candidatesN and
   --no-unedited is set.
3. Position ranges are 1-indexed inclusive.
4. Reference coordinates in modprofile match the BAM reference coordinates.
5. User input T/U is treated equivalently for profile comparisons.

Dependencies
------------
- pysam
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Set, Tuple

import pysam


class Region(NamedTuple):
    start: int
    end: int
    min_rate: Optional[float]        # None when candidates are user-supplied
    candidates: Optional[List[int]]  # None when using profile + min_rate; sorted, deduped


# ── CLI ───────────────────────────────────────────────────────────────────────


def normalize_base(base: str) -> str:
    """Normalize nucleotide input so T and U are treated equivalently."""
    b = base.strip().upper()
    if b == "U":
        return "T"
    return b


def strip_all_suffixes(path: str) -> str:
    """Return the filename stem with all dot-suffixes removed.

    e.g. sample.unique.sorted.bam -> sample
    """
    p = Path(path)
    while p.suffix:
        p = p.with_suffix("")
    return p.name


def detect_regions(argv: List[str]) -> int:
    """Scan argv for --startN / --endN / --min-rateN / --candidatesN and return the highest N found."""
    pattern = re.compile(r"^--(start|end|min-rate|candidates)(\d+)$")
    max_n = 0
    for arg in argv:
        m = pattern.match(arg)
        if m:
            max_n = max(max_n, int(m.group(2)))
    return max_n


def parse_args(argv: List[str]) -> argparse.Namespace:
    """Build and run the argument parser; requires at least one region to be present in argv."""
    max_n = detect_regions(argv)
    if max_n < 1:
        print(
            "Error: at least one region must be specified via "
            "--start1 --end1 --min-rate1  (or --candidates1 instead of --min-rate1)",
            file=sys.stderr,
        )
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description=(
            "Select QNAMEs from modified and untreated BAMs across one or more "
            "genomic regions. Threshold determined jointly from both BAMs per region."
        )
    )

    parser.add_argument(
        "--modbam", required=True,
        help="Modified condition BAM (e.g. ADAR-SHAPE). Coordinate-sorted and indexed."
    )
    parser.add_argument(
        "--untbam", required=True,
        help="Untreated condition BAM (e.g. ADAR-DMSO). Coordinate-sorted and indexed."
    )
    parser.add_argument(
        "--modprofile", required=False, default=None,
        help=(
            "Profile.txt from the modified condition (tab-delimited). Required unless "
            "all regions use --candidatesN and --no-unedited is set."
        )
    )
    parser.add_argument(
        "--refbase", default="A",
        help="Reference base to select from profile Sequence column. Default: A"
    )
    parser.add_argument(
        "--readbase", default="G",
        help="Read base required in alignments at selected positions. Default: G"
    )
    parser.add_argument(
        "--edited-min-depth", type=int, default=10000,
        help="Minimum number of QNAMEs required in both BAMs. Default: 10000"
    )
    parser.add_argument(
        "--minimum-coverage", type=float, default=75,
        help=(
            "Minimum %% of region length each read's aligned span must cover "
            "to be included in any counting step. Default: 75"
        )
    )
    parser.add_argument(
        "--contig", default=None,
        help=(
            "Reference/contig name to fetch from the BAMs. "
            "Auto-inferred if exactly one contig exists."
        )
    )
    parser.add_argument(
        "--out-dir", default=None,
        help="Directory for output files. Default: current working directory."
    )
    parser.add_argument(
        "--no-unedited", action="store_true", default=False,
        help="Skip generation of unmodified QNAME lists. Default: unmodified lists are produced."
    )
    parser.add_argument(
        "--unedited-min-depth", type=int, default=10000,
        help="Minimum reads required in both BAMs for unmodified threshold selection. Default: 10000"
    )
    parser.add_argument(
        "--max-tolerance", type=int, default=3,
        help=(
            "Maximum non-candidate A→G mutations allowed per molecule. "
            "Tolerance levels 0..N are each evaluated in a single BAM pass; "
            "the lowest passing level is selected. Default: 3"
        )
    )

    for n in range(1, max_n + 1):
        parser.add_argument(f"--start{n}", type=int, default=None,
                            help=f"1-indexed inclusive start for region {n}.")
        parser.add_argument(f"--end{n}", type=int, default=None,
                            help=f"1-indexed inclusive end for region {n}.")
        parser.add_argument(f"--min-rate{n}", type=float, default=None,
                            help=f"Min Modified_rate threshold for region {n}.")
        parser.add_argument(f"--candidates{n}", type=str, default=None,
                            help=(
                                f"Comma-separated 1-indexed candidate positions for region {n}. "
                                f"Mutually exclusive with --min-rate{n}."
                            ))

    return parser.parse_args(argv)


def parse_candidate_positions_str(s: str, start: int, end: int, region_n: int) -> List[int]:
    """Parse, deduplicate, sort, and range-validate a comma-separated position string."""
    try:
        positions = [int(x.strip()) for x in s.split(",") if x.strip()]
    except ValueError as exc:
        raise ValueError(
            f"Region {region_n}: --candidates{region_n} contains a non-integer value: {exc}"
        )
    if not positions:
        raise ValueError(f"Region {region_n}: --candidates{region_n} produced no positions")
    out_of_range = [p for p in positions if p < start or p > end]
    if out_of_range:
        raise ValueError(
            f"Region {region_n}: positions {out_of_range} in --candidates{region_n} "
            f"fall outside region [{start}, {end}]"
        )
    return sorted(set(positions))


def collect_regions(args: argparse.Namespace, max_n: int) -> List[Region]:
    """Build Region objects from parsed args; enforces mutual exclusivity of --min-rateN and --candidatesN."""
    regions = []
    for n in range(1, max_n + 1):
        start = getattr(args, f"start{n}", None)
        end = getattr(args, f"end{n}", None)
        min_rate = getattr(args, f"min_rate{n}", None)
        candidates_str = getattr(args, f"candidates{n}", None)

        location_provided = [x is not None for x in (start, end)]
        if not any(location_provided):
            continue  # slot N unused entirely
        if not all(location_provided):
            missing = [
                f"--{name}{n}"
                for name, v in zip(("start", "end"), (start, end))
                if v is None
            ]
            raise ValueError(f"Region {n} is incomplete — missing: {', '.join(missing)}")

        if min_rate is not None and candidates_str is not None:
            raise ValueError(
                f"Region {n}: --min-rate{n} and --candidates{n} are mutually exclusive"
            )
        if min_rate is None and candidates_str is None:
            raise ValueError(
                f"Region {n}: must provide either --min-rate{n} or --candidates{n}"
            )

        if candidates_str is not None:
            candidates = parse_candidate_positions_str(candidates_str, start, end, n)
            regions.append(Region(start=start, end=end, min_rate=None, candidates=candidates))
        else:
            regions.append(Region(start=start, end=end, min_rate=min_rate, candidates=None))

    if not regions:
        raise ValueError("No complete regions found.")
    return regions


def validate_args(args: argparse.Namespace, regions: List[Region]) -> None:
    """Validate all CLI arguments; normalizes refbase/readbase in-place on args."""
    args.refbase = normalize_base(args.refbase)
    args.readbase = normalize_base(args.readbase)

    if args.refbase not in {"A", "C", "G", "T"}:
        raise ValueError("--refbase must be one of A, C, G, T, U")
    if args.readbase not in {"A", "C", "G", "T"}:
        raise ValueError("--readbase must be one of A, C, G, T, U")
    if args.edited_min_depth < 1:
        raise ValueError("--edited-min-depth must be >= 1")
    if not (0 < args.minimum_coverage <= 100):
        raise ValueError("--minimum-coverage must be > 0 and <= 100")
    if args.unedited_min_depth < 1:
        raise ValueError("--unedited-min-depth must be >= 1")
    if args.max_tolerance < 0:
        raise ValueError("--max-tolerance must be >= 0")

    for bam_arg, label in [(args.modbam, "--modbam"), (args.untbam, "--untbam")]:
        if not Path(bam_arg).is_file():
            raise FileNotFoundError(f"{label} not found: {bam_arg}")
        if not Path(f"{bam_arg}.bai").is_file():
            raise FileNotFoundError(f"{label} index not found: {bam_arg}.bai")

    needs_profile = (not args.no_unedited) or any(r.min_rate is not None for r in regions)
    if needs_profile:
        if args.modprofile is None:
            raise ValueError(
                "--modprofile is required unless all regions use --candidatesN "
                "and --no-unedited is set"
            )
        if not Path(args.modprofile).is_file():
            raise FileNotFoundError(f"--modprofile not found: {args.modprofile}")

    if args.out_dir is not None and not Path(args.out_dir).is_dir():
        raise FileNotFoundError(f"--out-dir not found: {args.out_dir}")

    for i, r in enumerate(regions, 1):
        if r.start < 1:
            raise ValueError(f"Region {i}: --start{i} must be >= 1")
        if r.end < r.start:
            raise ValueError(f"Region {i}: --end{i} must be >= --start{i}")
        if r.min_rate is not None and r.min_rate <= 0:
            raise ValueError(f"Region {i}: --min-rate{i} must be > 0")


# ── Profile reading ───────────────────────────────────────────────────────────


def get_candidate_positions(
    profile_path: str,
    start: int,
    end: int,
    refbase: str,
    min_rate: float,
) -> List[int]:
    """Return 1-indexed positions in [start, end] where sequence == refbase and Modified_rate >= min_rate."""
    positions: List[int] = []

    with open(profile_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("profile.txt appears to have no header row")

        required = {"Nucleotide", "Sequence", "Modified_rate"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"profile.txt is missing required columns: {sorted(missing)}")

        for row in reader:
            try:
                pos = int(row["Nucleotide"])
            except Exception as exc:
                raise ValueError(
                    f"Could not parse Nucleotide value: {row.get('Nucleotide')!r}"
                ) from exc

            if pos < start or pos > end:
                continue

            seq_base = normalize_base(row["Sequence"])
            if seq_base != refbase:
                continue

            try:
                mod_rate = float(row["Modified_rate"])
            except Exception as exc:
                raise ValueError(
                    f"Could not parse Modified_rate value: {row.get('Modified_rate')!r}"
                ) from exc

            if mod_rate >= min_rate:
                positions.append(pos)

    return positions


def get_unmodified_positions(
    profile_path: str,
    start: int,
    end: int,
    refbase: str,
) -> List[int]:
    """Return all 1-indexed positions in [start, end] where sequence == refbase, regardless of Modified_rate."""
    positions: List[int] = []

    with open(profile_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("profile.txt appears to have no header row")

        required = {"Nucleotide", "Sequence"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"profile.txt is missing required columns: {sorted(missing)}")

        for row in reader:
            try:
                pos = int(row["Nucleotide"])
            except Exception as exc:
                raise ValueError(
                    f"Could not parse Nucleotide value: {row.get('Nucleotide')!r}"
                ) from exc

            if pos < start or pos > end:
                continue

            if normalize_base(row["Sequence"]) == refbase:
                positions.append(pos)

    return positions


# ── BAM scanning ──────────────────────────────────────────────────────────────


def is_primary_mapped(read: pysam.AlignedSegment) -> bool:
    """Return True if the read is mapped, non-secondary, and non-supplementary."""
    return (
        (not read.is_unmapped)
        and (not read.is_secondary)
        and (not read.is_supplementary)
    )


def union_coverage(intervals: List[Tuple[int, int]]) -> int:
    """Return the total non-overlapping length of a list of [start, end) intervals."""
    if not intervals:
        return 0
    total = 0
    merged = sorted(intervals)
    cur_start, cur_end = merged[0]
    for s, e in merged[1:]:
        if s <= cur_end:
            cur_end = max(cur_end, e)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = s, e
    total += cur_end - cur_start
    return total


def infer_contig_name(bam_path: str, user_contig: Optional[str]) -> str:
    """Return the contig to use: user-supplied if given, or the sole contig if the BAM has exactly one."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = list(bam.references)

    if user_contig is not None:
        if user_contig not in refs:
            raise ValueError(
                f"--contig {user_contig!r} not found in BAM header. "
                f"Available contigs: {refs}"
            )
        return user_contig

    if len(refs) == 1:
        return refs[0]

    raise ValueError(
        "BAM contains multiple contigs. Please provide --contig explicitly. "
        f"Available contigs: {refs}"
    )


def build_qname_to_matchcount(
    bam_path: str,
    contig: str,
    start: int,
    end: int,
    candidate_positions: Sequence[int],
    readbase: str,
    minimum_coverage: float,
) -> Dict[str, int]:
    """
    Scan the BAM and return {QNAME: n_distinct_candidate_positions_showing_readbase}
    for all coverage-qualified molecules. Uses union-of-mates position and interval
    tracking so concordant pairs are handled correctly.
    """
    candidate_set = set(candidate_positions)
    fetch_start = start - 1
    fetch_end = end
    min_span = minimum_coverage / 100.0 * (end - start + 1)
    qname_to_intervals: Dict[str, List[Tuple[int, int]]] = {}
    qname_to_hit_positions: Dict[str, Set[int]] = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(contig, fetch_start, fetch_end):
            if not is_primary_mapped(read):
                continue
            if read.query_sequence is None:
                continue

            qname = read.query_name

            mate_start = max(read.reference_start, fetch_start)
            mate_end   = min(read.reference_end,   fetch_end)
            if mate_end > mate_start:
                if qname not in qname_to_intervals:
                    qname_to_intervals[qname] = []
                qname_to_intervals[qname].append((mate_start, mate_end))

            if qname not in qname_to_hit_positions:
                qname_to_hit_positions[qname] = set()

            seq = read.query_sequence.upper()
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                if query_pos is None or ref_pos is None:
                    continue
                ref_pos_1 = ref_pos + 1
                if ref_pos_1 in candidate_set and seq[query_pos] == readbase:
                    qname_to_hit_positions[qname].add(ref_pos_1)

    coverage_qualified = {
        q for q, ivs in qname_to_intervals.items()
        if union_coverage(ivs) >= min_span
    }
    return {
        q: len(positions)
        for q, positions in qname_to_hit_positions.items()
        if q in coverage_qualified
    }


def build_unmodified_molecule_data(
    bam_path: str,
    contig: str,
    start: int,
    end: int,
    unmodified_positions: Sequence[int],
    candidate_set: set,
    refbase: str,
    readbase: str,
    minimum_coverage: float,
) -> Dict[str, Tuple[int, int]]:
    """
    Single BAM pass. Returns {QNAME: (n_noncand_muts, n_refbase_covered)} for
    coverage-qualified molecules with no G at any candidate position.

    n_noncand_muts   = count of distinct non-candidate positions showing readbase
                       (union across all mates).
    n_refbase_covered = count of distinct positions showing refbase (union across
                       all mates).

    Tolerance filtering is deferred to post-scan (see derive_per_tol_matchcounts).
    """
    unmod_set = set(unmodified_positions)
    fetch_start = start - 1
    fetch_end = end
    min_span = minimum_coverage / 100.0 * (end - start + 1)

    qname_to_intervals: Dict[str, List[Tuple[int, int]]] = {}
    cand_disqualified: Set[str] = set()
    qname_to_noncand_mut_pos: Dict[str, Set[int]] = {}
    qname_to_refbase_pos: Dict[str, Set[int]] = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(contig, fetch_start, fetch_end):
            if not is_primary_mapped(read):
                continue
            if read.query_sequence is None:
                continue

            qname = read.query_name

            mate_start = max(read.reference_start, fetch_start)
            mate_end   = min(read.reference_end,   fetch_end)
            if mate_end > mate_start:
                if qname not in qname_to_intervals:
                    qname_to_intervals[qname] = []
                qname_to_intervals[qname].append((mate_start, mate_end))

            if qname in cand_disqualified:
                continue

            seq = read.query_sequence.upper()
            covered: Dict[int, str] = {}
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                if query_pos is None or ref_pos is None:
                    continue
                ref_pos_1 = ref_pos + 1
                if ref_pos_1 in unmod_set:
                    covered[ref_pos_1] = seq[query_pos]

            if not covered:
                continue

            # Any readbase at a candidate position disqualifies the entire molecule
            if any(base == readbase and pos in candidate_set
                   for pos, base in covered.items()):
                cand_disqualified.add(qname)
                continue

            # Initialise accumulators on first accepted mate for this QNAME
            if qname not in qname_to_noncand_mut_pos:
                qname_to_noncand_mut_pos[qname] = set()
                qname_to_refbase_pos[qname] = set()

            # Accumulate non-candidate mutation positions (union across mates)
            for pos, base in covered.items():
                if base == readbase and pos not in candidate_set:
                    qname_to_noncand_mut_pos[qname].add(pos)

            # Accumulate refbase coverage positions (union across mates)
            for pos, base in covered.items():
                if base == refbase:
                    qname_to_refbase_pos[qname].add(pos)

    coverage_qualified = {
        q for q, ivs in qname_to_intervals.items()
        if union_coverage(ivs) >= min_span
    }
    return {
        qname: (len(qname_to_noncand_mut_pos[qname]), len(positions))
        for qname, positions in qname_to_refbase_pos.items()
        if qname not in cand_disqualified and qname in coverage_qualified
    }


def derive_per_tol_matchcounts(
    molecule_data: Dict[str, Tuple[int, int]],
    max_tolerance: int,
) -> Dict[int, Dict[str, int]]:
    """
    Returns {tol: {QNAME: n_refbase_covered}} for tol in 0..max_tolerance.
    A molecule qualifies at tol t if its n_noncand_muts <= t.
    """
    return {
        t: {q: mc for q, (nm, mc) in molecule_data.items() if nm <= t}
        for t in range(max_tolerance + 1)
    }


# ── Threshold selection & output ──────────────────────────────────────────────


def count_ge(matchcounts: Iterable[int], threshold: int) -> int:
    """Count how many values in matchcounts are >= threshold."""
    return sum(1 for x in matchcounts if x >= threshold)


def choose_best_threshold(
    unt_matchcounts: Dict[str, int],
    mod_matchcounts: Dict[str, int],
    n_candidates: int,
    min_depth: int,
) -> int:
    """
    Choose the highest threshold T where both BAMs have >= min_depth reads with
    matchcount >= T. Falls back toward 1 if the initial threshold cannot be met.
    """
    if n_candidates < 1:
        raise ValueError("No candidate positions available; cannot choose threshold")

    unt_counts = list(unt_matchcounts.values())
    mod_counts = list(mod_matchcounts.values())

    if not unt_counts or not mod_counts:
        raise ValueError("No reads overlapped the requested region in one or both BAMs")

    def both_meet(k: int) -> bool:
        return count_ge(unt_counts, k) >= min_depth and count_ge(mod_counts, k) >= min_depth

    initial_threshold = n_candidates if n_candidates <= 10 else 10

    if both_meet(initial_threshold):
        best = initial_threshold
        for k in range(initial_threshold + 1, n_candidates + 1):
            if both_meet(k):
                best = k
            else:
                break
        return best

    for k in range(initial_threshold - 1, 0, -1):
        if both_meet(k):
            return k

    # Even the loosest threshold (>=1 edited position) does not reach min_depth reads
    # in both BAMs: the edited set will be sparse and the downstream re-ShapeMapper
    # reactivity profile largely undefined (low-depth positions are masked out).
    sys.stderr.write(
        f"WARNING: fewer than {min_depth} reads meet even the >=1 edited-position "
        f"threshold in both BAMs; the edited read set is very small. Consider a lower "
        f"--min-rate (or different --candidates), or deeper sequencing.\n"
    )
    return 1


def choose_unmod_tolerance(
    unt_per_tol: Dict[int, Dict[str, int]],
    mod_per_tol: Dict[int, Dict[str, int]],
    n_unmod: int,
    max_tolerance: int,
    min_depth: int,
    minimum_coverage: float,
) -> Tuple[int, int]:
    """
    Return (chosen_tol, threshold) where chosen_tol is the lowest t in 0..max_tolerance
    where both BAMs have >= min_depth reads with matchcount >= threshold.

    threshold = max(1, ceil(minimum_coverage / 100 * n_unmod))  — fixed across all t.

    Exit 1 with a diagnostic table if no tolerance level meets the criterion.
    """
    threshold = max(1, math.ceil(minimum_coverage / 100.0 * n_unmod))
    results: List[Tuple[int, int, int]] = []
    for t in range(max_tolerance + 1):
        unt_n = count_ge(unt_per_tol[t].values(), threshold)
        mod_n = count_ge(mod_per_tol[t].values(), threshold)
        results.append((t, unt_n, mod_n))
        if unt_n >= min_depth and mod_n >= min_depth:
            return t, threshold

    lines = [
        f"ERROR: no tolerance level 0..{max_tolerance} achieved {min_depth} reads "
        f"at threshold >= {threshold} of {n_unmod} positions "
        f"({minimum_coverage:.4g}% coverage).",
        f"  {'tol':>4}  {'unt_reads':>9}  {'mod_reads':>9}",
    ]
    for t, un, mn in results:
        lines.append(f"  {t:>4}  {un:>9}  {mn:>9}")
    lines.append("Increase --max-tolerance or reduce --unedited-min-depth.")
    sys.exit("\n".join(lines))


def get_final_qnames(
    qname_to_matchcount: Dict[str, int],
    threshold: int,
) -> List[str]:
    """Return a sorted list of QNAMEs whose matchcount meets or exceeds threshold."""
    return sorted([q for q, c in qname_to_matchcount.items() if c >= threshold])


def minrate_tag(min_rate: float) -> str:
    """Convert min_rate to a filename-safe tag with minimum 2 decimal digits."""
    s = str(min_rate)
    decimal = s.split(".")[-1] if "." in s else "0"
    if len(decimal) < 2:
        decimal = decimal.ljust(2, "0")
    return f"minrate{decimal}"


def bam_output_path(
    bam_stem: str,
    out_dir: Path,
    refbase: str,
    readbase: str,
    start: int,
    end: int,
    min_rate: Optional[float],
    threshold: int,
    n_candidates: int,
    min_depth: int,
) -> Path:
    """Return the output Path for a sorted-reads BAM (does not write)."""
    rate_tag = minrate_tag(min_rate) if min_rate is not None else "usercandidates"
    kreads = min_depth // 1000
    filename = (
        f"{bam_stem}_{refbase}to{readbase}_range{start}-{end}"
        f"_{rate_tag}_atleast{threshold}of{n_candidates}pos"
        f"_min{kreads}kreads.bam"
    )
    return out_dir / filename


def unmodified_bam_output_path(
    bam_stem: str,
    out_dir: Path,
    refbase: str,
    readbase: str,
    start: int,
    end: int,
    threshold: int,
    n_unmod: int,
    chosen_tol: int,
) -> Path:
    """Return the output Path for an unmodified-reads BAM (does not write)."""
    filename = (
        f"{bam_stem}_{refbase}to{readbase}_range{start}-{end}"
        f"_unmodified_atleast{threshold}of{n_unmod}pos"
        f"_tol{chosen_tol}.bam"
    )
    return out_dir / filename


def write_all_bams_single_pass(
    source_bam_path: str,
    outputs: List[Tuple[Set[str], Path]],
) -> None:
    """
    Write multiple output BAMs from a single sequential scan of the source BAM.
    More efficient than one scan per output file when multiple regions or output
    types are requested.
    """
    with pysam.AlignmentFile(source_bam_path, "rb") as src:
        writers = [pysam.AlignmentFile(str(p), "wb", template=src) for _, p in outputs]
        try:
            for read in src:
                qname = read.query_name
                for writer, (qname_set, _) in zip(writers, outputs):
                    if qname in qname_set:
                        writer.write(read)
        finally:
            for w in writers:
                w.close()
    for _, out_path in outputs:
        pysam.index(str(out_path))


# ── Logging ───────────────────────────────────────────────────────────────────


def region_mode_label(region: Region) -> str:
    """Return a short human-readable string describing the candidate-selection mode for a region."""
    if region.min_rate is not None:
        return f"min-rate={region.min_rate}"
    return "user-supplied candidates"


SEP = "=" * 80
SUBSEP = "  " + "-" * 76


def hist_table(matchcounts: Dict[str, int], n_candidates: int,
               top_n: Optional[int] = None) -> str:
    """
    Build a right-aligned 3-column table:
      # Positions | # Reads (exact) | # Reads (>= n)
    Covers all values from n_candidates down to 0.
    If top_n is set, only the top_n rows (highest position counts) are shown,
    followed by a note on the number of omitted rows.
    """
    exact = Counter(matchcounts.values())
    rows: List[tuple] = []
    cumulative = 0
    for k in range(n_candidates, -1, -1):
        e = exact.get(k, 0)
        cumulative += e
        rows.append((k, e, cumulative))

    total_rows = len(rows)
    if top_n is not None and top_n < total_rows:
        shown = rows[:top_n]
        omitted = total_rows - top_n
    else:
        shown = rows
        omitted = 0

    h1, h2, h3 = "# Positions", "# Reads (exact)", "# Reads (>= n)"
    w1 = max(len(h1), len(str(n_candidates)))
    w2 = max(len(h2), len(str(max(exact.values(), default=0))))
    w3 = max(len(h3), len(str(sum(exact.values()))))

    lines = [
        f"  {h1:>{w1}}   {h2:>{w2}}   {h3:>{w3}}",
        f"  {'-'*w1}   {'-'*w2}   {'-'*w3}",
    ]
    for k, e, c in shown:
        lines.append(f"  {k:>{w1}}   {e:>{w2}}   {c:>{w3}}")
    if omitted:
        lines.append(f"  ... ({omitted} rows with lower position counts omitted)")
    return "\n".join(lines)


def hist_table_multi_tol(
    per_tol: Dict[int, Dict[str, int]],
    n_unmod: int,
    max_tolerance: int,
    threshold: int,
) -> str:
    """
    Multi-tolerance histogram. Rows: n_unmod down to 0 (no cap).
    Columns: # Positions, then for each tol t: (tol{t} exact, tol{t} >=n).
    The row at k == threshold is marked with  *  to indicate the selection cutoff.
    """
    # Pre-compute exact counts and cumulative sums per tol
    per_exact: Dict[int, Counter] = {}
    per_cum: Dict[int, Dict[int, int]] = {}
    for t in range(max_tolerance + 1):
        exact = Counter(per_tol[t].values())
        per_exact[t] = exact
        cum: Dict[int, int] = {}
        running = 0
        for k in range(n_unmod, -1, -1):
            running += exact.get(k, 0)
            cum[k] = running
        per_cum[t] = cum

    max_count = max(
        (sum(per_exact[t].values()) for t in range(max_tolerance + 1)),
        default=0,
    )

    h_pos = "# Positions"
    w_pos = max(len(h_pos), len(str(n_unmod)))
    w_e = max(len(f"tol{max_tolerance} exact"), len(str(max_count)))
    w_c = max(len(f"tol{max_tolerance} >=n"), len(str(max_count)))
    w_mark = 2  # space for "  *" or "   "

    # Header
    header = f"  {h_pos:>{w_pos}}"
    for t in range(max_tolerance + 1):
        header += f"   {f'tol{t} exact':>{w_e}}   {f'tol{t} >=n':>{w_c}}"
    sep = f"  {'-'*w_pos}"
    for _ in range(max_tolerance + 1):
        sep += f"   {'-'*w_e}   {'-'*w_c}"

    lines = [header, sep]
    for k in range(n_unmod, -1, -1):
        row = f"  {k:>{w_pos}}"
        for t in range(max_tolerance + 1):
            row += f"   {per_exact[t].get(k, 0):>{w_e}}   {per_cum[t].get(k, 0):>{w_c}}"
        row += "  *" if k == threshold else ""
        lines.append(row)

    return "\n".join(lines)


def log_bam_section(log, label: str, out_path: str, qname_count: int,
                    matchcounts: Dict[str, int], n_candidates: int) -> None:
    """Write one BAM section (output path, QNAME count, edited-reads histogram) to log."""
    log.write(f"{SUBSEP}\n")
    log.write(f"  {label}\n")
    log.write(f"{SUBSEP}\n")
    log.write(f"  Output:  {out_path}\n")
    log.write(f"  QNAMEs:  {qname_count}\n\n")
    log.write(hist_table(matchcounts, n_candidates))
    log.write("\n\n")


def log_unmod_section(
    log,
    refbase: str,
    n_unmod: int,
    max_tolerance: int,
    threshold: int,
    chosen_tol: int,
    unt_per_tol: Dict[int, Dict[str, int]],
    mod_per_tol: Dict[int, Dict[str, int]],
    unt_unmod_out: str,
    mod_unmod_out: str,
    unt_unmod_count: int,
    mod_unmod_count: int,
    minimum_coverage: float,
) -> None:
    """Write the combined unmodified-reads section (both BAMs, multi-tolerance histogram)."""
    unt_chosen = count_ge(unt_per_tol[chosen_tol].values(), threshold)
    mod_chosen = count_ge(mod_per_tol[chosen_tol].values(), threshold)

    log.write(f"{SUBSEP}\n")
    log.write(f"  Unmodified reads  (max-tolerance={max_tolerance})\n")
    log.write(f"{SUBSEP}\n")
    log.write(f"  Unmodified positions (all {refbase} in region, no min-rate filter): {n_unmod}\n")
    log.write(
        f"  Chosen:  tol{chosen_tol}  "
        f"(threshold >= {threshold} of {n_unmod} positions "
        f"[{minimum_coverage:.4g}%],  unt={unt_chosen},  mod={mod_chosen})\n\n"
    )
    log.write("  Untreated BAM:\n")
    log.write(hist_table_multi_tol(unt_per_tol, n_unmod, max_tolerance, threshold))
    log.write("\n\n")
    log.write("  Modified BAM:\n")
    log.write(hist_table_multi_tol(mod_per_tol, n_unmod, max_tolerance, threshold))
    log.write("\n\n")
    log.write(f"  Untreated unmodified output:  {unt_unmod_out}  ({unt_unmod_count} QNAMEs)\n")
    log.write(f"  Modified  unmodified output:  {mod_unmod_out}  ({mod_unmod_count} QNAMEs)\n")
    log.write("\n")


def progress(msg: str) -> None:
    """Print a progress message to stderr."""
    print(msg, file=sys.stderr)


def log_run_parameters(log, args, contig, out_dir, regions) -> None:
    """Write the run-parameter header block to the log file."""
    log.write(f"{SEP}\n")
    log.write("Run Parameters\n")
    log.write(f"{SEP}\n")
    log.write(f"Date/time:              {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log.write(f"--modbam:               {args.modbam}\n")
    log.write(f"--untbam:               {args.untbam}\n")
    log.write(f"--modprofile:           {args.modprofile if args.modprofile else 'N/A'}\n")
    log.write(f"--refbase:              {args.refbase}\n")
    log.write(f"--readbase:             {args.readbase}\n")
    log.write(f"--edited-min-depth:            {args.edited_min_depth}\n")
    log.write(f"--minimum-coverage:     {args.minimum_coverage}\n")
    log.write(f"--contig:               {contig}\n")
    log.write(f"--out-dir:              {out_dir}\n\n")
    log.write("Regions:\n")
    for i, r in enumerate(regions, 1):
        log.write(f"  Region {i}: start={r.start}  end={r.end}  {region_mode_label(r)}\n")
    log.write(f"--no-unedited:        {args.no_unedited}\n")
    if not args.no_unedited:
        log.write(f"--unedited-min-depth: {args.unedited_min_depth}\n")
        log.write(f"--max-tolerance:        {args.max_tolerance}\n")
    log.write("\n")


def log_region_skipped(log, i: int, region: Region, msg: str) -> None:
    """Write a skipped-region entry to the log file."""
    log.write(f"{SEP}\n")
    log.write(f"Region {i}: {region.start}-{region.end}  ({region_mode_label(region)})\n")
    log.write(f"{SEP}\n")
    log.write(f"SKIPPED: {msg}\n\n")


def log_region_candidates(log, i: int, region: Region,
                           candidate_positions: List[int], threshold: int) -> None:
    """Write the candidate positions and chosen threshold for a region to the log file."""
    n = len(candidate_positions)
    pos_list = ",".join(map(str, candidate_positions))
    log.write(f"{SEP}\n")
    log.write(f"Region {i}: {region.start}-{region.end}  ({region_mode_label(region)})\n")
    log.write(f"{SEP}\n")
    log.write(f"Candidate positions ({n}):\n  {pos_list}\n\n")
    log.write(f"Chosen threshold:  >= {threshold} of {n}  (determined from both BAMs)\n\n")


def log_region_footer(log, elapsed: float) -> None:
    """Write the per-region elapsed time to the log."""
    log.write(f"Elapsed (this region):  {elapsed:.2f}s\n\n")


def log_overall_summary(log, elapsed: float) -> None:
    """Write the overall elapsed time footer to the log file."""
    log.write(f"{SEP}\n")
    log.write("Overall\n")
    log.write(f"{SEP}\n")
    log.write(f"Total elapsed:  {elapsed:.2f}s\n")


# ── Entry point ───────────────────────────────────────────────────────────────


def main() -> None:
    overall_start = time.time()
    argv = sys.argv[1:]
    max_n = detect_regions(argv)

    args = parse_args(argv)
    regions = collect_regions(args, max_n)
    validate_args(args, regions)

    out_dir = Path(args.out_dir) if args.out_dir else Path.cwd()
    unt_stem = strip_all_suffixes(args.untbam)
    mod_stem = strip_all_suffixes(args.modbam)
    log_path = out_dir / f"{mod_stem}_log.txt"

    progress("Determining contig name...")
    contig = infer_contig_name(args.modbam, args.contig)

    unt_pending: List[Tuple[Set[str], Path]] = []
    mod_pending: List[Tuple[Set[str], Path]] = []

    with open(log_path, "w") as log:
        log_run_parameters(log, args, contig, out_dir, regions)

        for i, region in enumerate(regions, 1):
            region_start_time = time.time()
            progress(
                f"\n[Region {i}/{len(regions)}] {region.start}-{region.end} "
                f"({region_mode_label(region)})"
            )

            if region.candidates is not None:
                candidate_positions = region.candidates
            else:
                progress("  Reading modprofile and selecting candidate positions...")
                candidate_positions = get_candidate_positions(
                    profile_path=args.modprofile,
                    start=region.start,
                    end=region.end,
                    refbase=args.refbase,
                    min_rate=region.min_rate,
                )

            if not candidate_positions:
                msg = (
                    f"no candidate positions met the filters "
                    f"(range={region.start}-{region.end}, refbase={args.refbase}, "
                    f"{region_mode_label(region)})"
                )
                progress(f"  WARNING: {msg}")
                log_region_skipped(log, i, region, msg)
                continue

            progress("  Scanning untreated BAM...")
            unt_matchcounts = build_qname_to_matchcount(
                bam_path=args.untbam,
                contig=contig,
                start=region.start,
                end=region.end,
                candidate_positions=candidate_positions,
                readbase=args.readbase,
                minimum_coverage=args.minimum_coverage,
            )

            progress("  Scanning modified BAM...")
            mod_matchcounts = build_qname_to_matchcount(
                bam_path=args.modbam,
                contig=contig,
                start=region.start,
                end=region.end,
                candidate_positions=candidate_positions,
                readbase=args.readbase,
                minimum_coverage=args.minimum_coverage,
            )

            progress("  Choosing threshold from both BAMs...")
            threshold = choose_best_threshold(
                unt_matchcounts=unt_matchcounts,
                mod_matchcounts=mod_matchcounts,
                n_candidates=len(candidate_positions),
                min_depth=args.edited_min_depth,
            )

            unt_qnames = get_final_qnames(unt_matchcounts, threshold)
            mod_qnames = get_final_qnames(mod_matchcounts, threshold)

            shared_path_kwargs = dict(
                out_dir=out_dir,
                refbase=args.refbase,
                readbase=args.readbase,
                start=region.start,
                end=region.end,
                min_rate=region.min_rate,
                threshold=threshold,
                n_candidates=len(candidate_positions),
                min_depth=args.edited_min_depth,
            )
            unt_path = bam_output_path(bam_stem=unt_stem, **shared_path_kwargs)
            mod_path = bam_output_path(bam_stem=mod_stem, **shared_path_kwargs)
            unt_pending.append((set(unt_qnames), unt_path))
            mod_pending.append((set(mod_qnames), mod_path))
            unt_out = str(unt_path)
            mod_out = str(mod_path)

            n = len(candidate_positions)
            log_region_candidates(log, i, region, candidate_positions, threshold)

            if not args.no_unedited:
                progress("  Getting unmodified positions...")
                unmod_positions = get_unmodified_positions(
                    profile_path=args.modprofile,
                    start=region.start,
                    end=region.end,
                    refbase=args.refbase,
                )
                candidate_set = set(candidate_positions)
                n_unmod = len(unmod_positions)

                progress("  Scanning untreated BAM for unmodified reads (single pass)...")
                unt_mol_data = build_unmodified_molecule_data(
                    bam_path=args.untbam, contig=contig,
                    start=region.start, end=region.end,
                    unmodified_positions=unmod_positions,
                    candidate_set=candidate_set,
                    refbase=args.refbase, readbase=args.readbase,
                    minimum_coverage=args.minimum_coverage,
                )

                progress("  Scanning modified BAM for unmodified reads (single pass)...")
                mod_mol_data = build_unmodified_molecule_data(
                    bam_path=args.modbam, contig=contig,
                    start=region.start, end=region.end,
                    unmodified_positions=unmod_positions,
                    candidate_set=candidate_set,
                    refbase=args.refbase, readbase=args.readbase,
                    minimum_coverage=args.minimum_coverage,
                )

                unt_per_tol = derive_per_tol_matchcounts(unt_mol_data, args.max_tolerance)
                mod_per_tol = derive_per_tol_matchcounts(mod_mol_data, args.max_tolerance)

                progress("  Selecting tolerance level...")
                chosen_tol, unmod_threshold = choose_unmod_tolerance(
                    unt_per_tol=unt_per_tol,
                    mod_per_tol=mod_per_tol,
                    n_unmod=n_unmod,
                    max_tolerance=args.max_tolerance,
                    min_depth=args.unedited_min_depth,
                    minimum_coverage=args.minimum_coverage,
                )

                unt_unmod_qnames = get_final_qnames(unt_per_tol[chosen_tol], unmod_threshold)
                mod_unmod_qnames = get_final_qnames(mod_per_tol[chosen_tol], unmod_threshold)

                unmod_path_kwargs = dict(
                    out_dir=out_dir,
                    refbase=args.refbase, readbase=args.readbase,
                    start=region.start, end=region.end,
                    threshold=unmod_threshold, n_unmod=n_unmod,
                    chosen_tol=chosen_tol,
                )
                unt_unmod_path = unmodified_bam_output_path(bam_stem=unt_stem, **unmod_path_kwargs)
                mod_unmod_path = unmodified_bam_output_path(bam_stem=mod_stem, **unmod_path_kwargs)
                unt_pending.append((set(unt_unmod_qnames), unt_unmod_path))
                mod_pending.append((set(mod_unmod_qnames), mod_unmod_path))
                unt_unmod_out = str(unt_unmod_path)
                mod_unmod_out = str(mod_unmod_path)
            else:
                n_unmod = 0
                chosen_tol = None
                unmod_threshold = None
                unt_unmod_out = mod_unmod_out = None
                unt_unmod_qnames = mod_unmod_qnames = []
                unt_per_tol = mod_per_tol = None

            region_elapsed = time.time() - region_start_time

            log_bam_section(log, "Untreated BAM", unt_out, len(unt_qnames), unt_matchcounts, n)
            log_bam_section(log, "Modified BAM", mod_out, len(mod_qnames), mod_matchcounts, n)

            if not args.no_unedited:
                log_unmod_section(
                    log,
                    refbase=args.refbase,
                    n_unmod=n_unmod,
                    max_tolerance=args.max_tolerance,
                    threshold=unmod_threshold,
                    chosen_tol=chosen_tol,
                    unt_per_tol=unt_per_tol,
                    mod_per_tol=mod_per_tol,
                    unt_unmod_out=unt_unmod_out,
                    mod_unmod_out=mod_unmod_out,
                    unt_unmod_count=len(unt_unmod_qnames),
                    mod_unmod_count=len(mod_unmod_qnames),
                    minimum_coverage=args.minimum_coverage,
                )

            log_region_footer(log, region_elapsed)

        if unt_pending:
            progress("\nWriting untreated BAM outputs (single pass)...")
            write_all_bams_single_pass(args.untbam, unt_pending)
        if mod_pending:
            progress("Writing modified BAM outputs (single pass)...")
            write_all_bams_single_pass(args.modbam, mod_pending)

        log_overall_summary(log, time.time() - overall_start)

    progress(f"\nDone. Log written to: {log_path}")


if __name__ == "__main__":
    main()
