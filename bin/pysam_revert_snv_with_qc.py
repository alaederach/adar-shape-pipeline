#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set

import pysam


VALID_BASES = {"A", "C", "G", "T"}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Revert read bases back to the reference base for aligned mismatches "
            "matching --refbase in the reference and --readbase in the read. "
            "Requires MD tags in the input SAM/BAM."
        )
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input SAM or BAM file"
    )

    parser.add_argument(
        "--out", "--output", "-o",
        dest="output",
        help=(
            "Output SAM or BAM file. "
            "Default: [input_stem]_[refbase]to[readbase]_reverted.[same extension as input]"
        )
    )

    parser.add_argument(
        "--refbase",
        required=True,
        help="Reference base to match (A, C, G, or T)"
    )

    parser.add_argument(
        "--readbase",
        required=True,
        help="Read base to match and revert back to --refbase (A, C, G, or T)"
    )

    parser.add_argument(
        "--region-start",
        type=int,
        default=None,
        help="Sorting region start position (1-based inclusive). Used for QC in/out-of-region counts."
    )

    parser.add_argument(
        "--region-end",
        type=int,
        default=None,
        help="Sorting region end position (1-based inclusive). Used for QC in/out-of-region counts."
    )

    parser.add_argument(
        "--qc-out",
        default=None,
        help="Output path for per-QNAME QC TSV. If omitted, no QC file is written."
    )

    args = parser.parse_args()

    args.refbase = args.refbase.upper()
    args.readbase = args.readbase.upper()

    if args.refbase not in VALID_BASES:
        parser.error("--refbase must be one of: A, C, G, T")

    if args.readbase not in VALID_BASES:
        parser.error("--readbase must be one of: A, C, G, T")

    if args.refbase == args.readbase:
        parser.error("--refbase and --readbase must differ")

    return args


def infer_pysam_modes(infile: Path, outfile: Path):
    in_suffix = infile.suffix.lower()
    out_suffix = outfile.suffix.lower()

    if in_suffix == ".sam":
        in_mode = "r"
    elif in_suffix == ".bam":
        in_mode = "rb"
    else:
        raise ValueError(f"Unsupported input file extension: {infile}")

    if out_suffix == ".sam":
        out_mode = "w"
    elif out_suffix == ".bam":
        out_mode = "wb"
    else:
        raise ValueError(f"Unsupported output file extension: {outfile}")

    return in_mode, out_mode


def default_output_path(infile: Path, refbase: str, readbase: str) -> Path:
    suffix = infile.suffix.lower()

    if suffix not in {".sam", ".bam"}:
        raise ValueError(
            f"Unsupported input file extension for default output naming: {infile}"
        )

    tag = f"_{refbase}to{readbase}_reverted"
    return infile.with_name(f"{infile.stem}{tag}{suffix}")


def alignment_type(read) -> str:
    if read.is_secondary:
        return "secondary"
    if read.is_supplementary:
        return "supplementary"
    if read.is_unmapped:
        return "unmapped"
    return "primary"


def revert_read(read: pysam.AlignedSegment, refbase: str, readbase: str) -> List[int]:
    """Revert refbase→readbase mismatches in read in-place. Returns 0-based ref positions reverted."""
    aligned_pairs = read.get_aligned_pairs(with_seq=True)
    read_seq = list(read.query_sequence)
    read_qscores = read.query_qualities
    reverted = []

    for query_pos, ref_pos, ref_base in aligned_pairs:
        if query_pos is None or ref_pos is None or ref_base is None:
            continue
        if ref_base.upper() == refbase and read_seq[query_pos].upper() == readbase:
            read_seq[query_pos] = refbase
            reverted.append(ref_pos)

    if reverted:
        read.query_sequence = "".join(read_seq)
        read.query_qualities = read_qscores
        if read.has_tag("MD"):
            read.set_tag("MD", None)
        if read.has_tag("NM"):
            read.set_tag("NM", None)

    return reverted


def compute_qc_rows(
    qname_positions: Dict[str, Set[int]],
    region_0start: Optional[int],
    region_0end: Optional[int],
) -> List[Dict]:
    """Build one QC row per QNAME from the union of reverted positions across all mates.

    alignment_type and n_mates_reverted are left as None and must be filled in by the caller.
    """
    rows = []
    have_region = region_0start is not None and region_0end is not None
    for qname, positions in qname_positions.items():
        sorted_pos = sorted(positions)
        if have_region:
            n_in = sum(1 for p in sorted_pos if region_0start <= p < region_0end)
        else:
            n_in = 0
        n_out = len(sorted_pos) - n_in
        rows.append({
            "qname": qname,
            "alignment_type": None,
            "n_mates_reverted": None,
            "n_reverted_in_region": n_in,
            "n_reverted_out_of_region": n_out,
            "n_reverted_total": len(sorted_pos),
            "reverted_positions_0based": ",".join(str(p) for p in sorted_pos),
            "reverted_positions_1based": ",".join(str(p + 1) for p in sorted_pos),
        })
    return rows


def main():
    args = parse_args()

    infile = Path(args.input)
    outfile = Path(args.output) if args.output else default_output_path(
        infile, args.refbase, args.readbase
    )

    if not infile.is_file():
        sys.exit(f"ERROR: input file not found: {infile}")

    try:
        in_mode, out_mode = infer_pysam_modes(infile, outfile)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")

    do_qc = args.qc_out is not None
    have_region = args.region_start is not None and args.region_end is not None
    # Convert 1-based inclusive region coords to 0-based half-open for comparison with pysam positions
    region_0start = (args.region_start - 1) if have_region else None
    region_0end = args.region_end if have_region else None

    n_total = 0
    n_modified_reads = 0
    n_modified_bases = 0

    # QC accumulators (keyed by QNAME)
    qname_to_positions: Dict[str, Set[int]] = {}  # union of reverted 0-based ref positions
    qname_mate_count:   Dict[str, int]      = {}  # number of mates that had ≥1 reversion
    qname_aln_types:    Dict[str, List[str]]= {}  # alignment_type of each reverted mate
    all_qnames:         Set[str]            = set()

    with pysam.AlignmentFile(str(infile), in_mode) as in_fh:
        with pysam.AlignmentFile(str(outfile), out_mode, template=in_fh) as out_fh:
            for read in in_fh.fetch(until_eof=True):
                n_total += 1

                if read.query_sequence is None:
                    out_fh.write(read)
                    continue

                try:
                    positions = revert_read(read, args.refbase, args.readbase)
                except ValueError as e:
                    sys.exit(
                        "ERROR: input appears to lack MD tags required for "
                        f"get_aligned_pairs(with_seq=True).\nUnderlying error: {e}"
                    )

                n_modified_bases += len(positions)
                if positions:
                    n_modified_reads += 1

                out_fh.write(read)

                if do_qc:
                    qname = read.query_name
                    all_qnames.add(qname)
                    if positions:
                        qname_to_positions.setdefault(qname, set()).update(positions)
                        qname_mate_count[qname] = qname_mate_count.get(qname, 0) + 1
                        qname_aln_types.setdefault(qname, []).append(alignment_type(read))

    if do_qc:
        qc_rows = compute_qc_rows(qname_to_positions, region_0start, region_0end)
        for row in qc_rows:
            qname = row["qname"]
            row["n_mates_reverted"] = qname_mate_count[qname]
            row["alignment_type"] = ",".join(qname_aln_types[qname])
        n_unmodified_qnames = len(all_qnames) - len(qname_to_positions)

        qc_path = Path(args.qc_out)
        qc_path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = [
            "qname",
            "alignment_type",
            "n_mates_reverted",
            "n_reverted_in_region",
            "n_reverted_out_of_region",
            "n_reverted_total",
            "reverted_positions_0based",
            "reverted_positions_1based",
        ]
        with open(qc_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(qc_rows)
        region_str = (
            f"{args.region_start}-{args.region_end} (1-based inclusive; "
            f"0-based half-open: [{region_0start}, {region_0end}))"
            if have_region else "not provided"
        )
        print(f"QC output: {qc_path}", file=sys.stderr)
        print(f"  Sorting region: {region_str}", file=sys.stderr)
        print(f"  QNAMEs with reversions (rows in QC file): {len(qc_rows):,}", file=sys.stderr)
        print(f"  QNAMEs with zero reversions (not in QC file): {n_unmodified_qnames:,}", file=sys.stderr)

    print(f"Input: {infile}", file=sys.stderr)
    print(f"Output: {outfile}", file=sys.stderr)
    print(f"Reference base matched: {args.refbase}", file=sys.stderr)
    print(f"Read base reverted: {args.readbase} -> {args.refbase}", file=sys.stderr)
    print(f"Reads processed: {n_total:,}", file=sys.stderr)
    print(f"Reads modified: {n_modified_reads:,}", file=sys.stderr)
    print(f"Bases reverted: {n_modified_bases:,}", file=sys.stderr)
    print()
    print()

if __name__ == "__main__":
    main()
