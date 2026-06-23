#!/usr/bin/env bash
# Fetch the bundled demo dataset (downsampled STK4 AluSq2 FASTQs) and unpack it
# into data/AluSq2/fastqs/. The FASTQs are hosted on UNC Dataverse rather than
# committed to git (they are large binaries; full-depth data goes to GEO/SRA on
# publication — see README).
#
# Usage:  bin/download_data.sh
set -euo pipefail

# UNC Dataverse — "download all files as zip" endpoint for the demo dataset.
#   DOI: https://doi.org/10.15139/S3/E2BTCM
# (override at runtime with ADAR_SHAPE_DATA_URL=... bin/download_data.sh)
DOI="${ADAR_SHAPE_DATA_DOI:-doi:10.15139/S3/E2BTCM}"
DATA_URL="${ADAR_SHAPE_DATA_URL:-https://dataverse.unc.edu/api/access/dataset/:persistentId?persistentId=${DOI}}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$ROOT/data"
ZIP="$DEST/adar-shape-demo-data.zip"

mkdir -p "$DEST"
echo "Downloading demo data from UNC Dataverse ($DOI) ..."
curl -fL --retry 3 -o "$ZIP" "$DATA_URL"

echo "Unpacking into $DEST/AluSq2 ..."
unzip -oq "$ZIP" -d "$DEST"
rm -f "$ZIP" "$DEST/MANIFEST.TXT"   # MANIFEST.TXT is added by Dataverse's bulk download

n=$(find "$DEST/AluSq2/fastqs" -name '*.fastq.gz' | wc -l | tr -d ' ')
echo "Done — $n FASTQ files in $DEST/AluSq2/fastqs/"
[ "$n" -eq 12 ] || { echo "WARNING: expected 12 FASTQs, found $n" >&2; exit 1; }
