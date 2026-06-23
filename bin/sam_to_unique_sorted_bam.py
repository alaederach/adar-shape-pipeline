#!/usr/bin/env python3

"""
Filter a SAM file to uniquely mapped reads and write a coordinate-sorted,
indexed BAM.

Output is written to the same directory as the input SAM, named:
    <stem>.unique.sorted.bam

Uniquely mapped is defined as:
  - Unpaired (merged) reads: FLAG 0x1 not set, no XS:i: tag
  - Concordant pairs: FLAG 0x1 set, FLAG 0x2 set, no XS:i: tag

The XS:i: tag is set by bowtie2 on any read with a suboptimal alignment,
indicating it could map to at least one other location. Reads without XS
mapped to exactly one location.

This approach replaces the name-sort + QNAME-grouping strategy used in v1,
allowing per-record filtering on a single pass through the SAM — no
name-sort step required.

Dependencies
------------
- pysam
- samtools
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Tuple

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a SAM to uniquely mapped reads and write a coordinate-sorted "
            "indexed BAM named <stem>.unique.sorted.bam."
        )
    )
    parser.add_argument("--sam", required=True, help="Input SAM file.")
    parser.add_argument(
        "--out",
        default=None,
        help=(
            "Output BAM path. Default: <sam_dir>/<stem>.unique.sorted.bam"
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads to pass to samtools sort. Default: 4",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    sam_path = Path(args.sam)
    if not sam_path.is_file():
        raise FileNotFoundError(f"Input SAM not found: {sam_path}")
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools not found in PATH")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")


def is_unique_mapper(read: pysam.AlignedSegment) -> bool:
    """
    Return True if the read is a uniquely mapped unpaired read or one mate
    of a concordantly aligned pair, with no alternative alignment location.

    Criteria:
      - Not unmapped, secondary, or supplementary
      - No XS:i: tag (no suboptimal alignment exists)
      - Either: unpaired (FLAG 0x1 not set)
      - Or:     concordant pair (FLAG 0x1 and FLAG 0x2 both set)
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False

    if read.has_tag("XS"):
        return False

    # Unpaired (merged) read
    if not read.is_paired:
        return True

    # Concordant pair mate
    if read.is_proper_pair:
        return True

    return False


def build_unique_sorted_bam(
    sam_path: str,
    out_bam: str,
    threads: int = 4,
) -> Tuple[int, int]:
    """
    Stream pipeline:
        pysam read SAM -> per-record filter -> samtools sort (coord) -> out_bam

    Returns (kept_records, removed_records).
    """
    coord_sort_cmd = [
        "samtools", "sort",
        "-O", "BAM",
        "-@", str(threads),
        "-o", out_bam,
        "-",
    ]

    p = subprocess.Popen(coord_sort_cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.stdin is None:
        raise RuntimeError("Failed to open stdin pipe to samtools coordinate sort")

    kept_records = 0
    removed_records = 0

    try:
        with pysam.AlignmentFile(sam_path, "r") as in_sam:
            with pysam.AlignmentFile(p.stdin, "wb", header=in_sam.header) as out_stream:
                for read in in_sam:
                    if is_unique_mapper(read):
                        out_stream.write(read)
                        kept_records += 1
                    else:
                        removed_records += 1
    finally:
        if p.stdin is not None:
            p.stdin.close()

    stderr = p.stderr.read().decode("utf-8", errors="replace") if p.stderr else ""
    rc = p.wait()

    if rc != 0:
        raise RuntimeError(f"samtools coordinate sort failed with exit code {rc}\n{stderr}")

    pysam.index(out_bam)

    return kept_records, removed_records


def main() -> None:
    start_time = time.time()
    args = parse_args()
    validate_args(args)

    sam_path = Path(args.sam)
    out_bam = args.out if args.out else str(sam_path.parent / f"{sam_path.stem}.unique.sorted.bam")

    print(f"Input SAM:  {sam_path}", file=sys.stderr)
    print(f"Output BAM: {out_bam}", file=sys.stderr)
    print("Filtering and coordinate-sorting...", file=sys.stderr)

    kept, removed = build_unique_sorted_bam(
        sam_path=str(sam_path),
        out_bam=out_bam,
        threads=args.threads,
    )

    elapsed = time.time() - start_time
    print("\nDone.")
    print(f"Output BAM:         {out_bam}")
    print(f"Kept alignments:    {kept}")
    print(f"Removed alignments: {removed}")
    print(f"Elapsed seconds:    {elapsed:.2f}")


if __name__ == "__main__":
    main()
