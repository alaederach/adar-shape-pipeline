#!/usr/bin/env bash
# Reproducibly build the bundled demo dataset (data/AluSq2/) by downsampling the
# full STK4 AluSq2 sequencing data to a uniform target read count per FASTQ.
#
# Each .fastq.gz is reduced to ~TARGET read pairs by keeping every Nth 4-line
# record (N = total/TARGET). The selection is positional, so paired R1/R2 files —
# which share the same record count — keep the SAME records and stay in sync.
# Whole records are preserved. FASTA + primers are copied unchanged.
#
# Usage:  bin/make_demo_data.sh <full_AluSq2_dir> [TARGET]
#   <full_AluSq2_dir>  dir containing: STK4-singleAlus-AluSq2-corrected.fa,
#                      primers/STK4_AluSq2_primers.txt,
#                      fastqs/{ADAR_SHAPE,ADAR_DMSO,MOCK_SHAPE,MOCK_DMSO}/rep1/*.fastq.gz
#   [TARGET]           target read pairs per file (default 100000)
set -euo pipefail

SRC="${1:?usage: make_demo_data.sh <full_AluSq2_dir> [TARGET]}"
TARGET="${2:-100000}"
DST="$(cd "$(dirname "$0")/.." && pwd)/data/AluSq2"

mkdir -p "$DST/primers"
cp "$SRC/STK4-singleAlus-AluSq2-corrected.fa" "$DST/"
cp "$SRC/primers/STK4_AluSq2_primers.txt"     "$DST/primers/"

find "$SRC/fastqs" -name '*.fastq.gz' | while read -r f; do
    rel="${f#"$SRC"/fastqs/}"
    out="$DST/fastqs/$rel"
    mkdir -p "$(dirname "$out")"
    total=$(( $(gzip -dc "$f" | wc -l) / 4 ))
    stride=$(( total / TARGET )); [ "$stride" -lt 1 ] && stride=1
    gzip -dc "$f" | awk -v s="$stride" 'int((NR-1)/4) % s == 0' | gzip > "$out"
    echo "  $rel : $total -> ~$(( total / stride )) reads (stride $stride)"
done
echo "demo data written to $DST"
