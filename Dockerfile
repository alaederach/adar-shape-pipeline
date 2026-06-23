# ShapeMapper2 (x86-64 Linux only) packaged so it can be called on non-Linux hosts.
# ShapeMapper2 has no macOS build; this image lets Nextflow run it on macOS — including
# Apple Silicon, under Docker Desktop's Rosetta-backed linux/amd64 emulation.
#
# Build:
#   docker build --platform linux/amd64 -t shape-adar/shapemapper:2.2.0 .
#
# It is invoked transparently by bin/shapemapper-docker, which params.shapemapper_bin
# points at. On a native Linux host you do not need this image at all — install
# ShapeMapper2 normally and set params.shapemapper_bin to the real `shapemapper`.

FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates bzip2 procps libgomp1 default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Pinned to v2.2 for consistency with the reference pipeline (the shapemapper.nf module's
# output naming, e.g. {name}_aligned.sam, matches 2.2; 2.3 renamed those outputs).
# The release tarball is packaged on macOS and may ship AppleDouble "._*" resource-fork
# files; Python's site.py would parse the bundled "._*.pth" ones as real .pth files and
# crash on their binary header (UnicodeDecodeError). Strip all "._*" files after extraction.
ARG SM_URL=https://github.com/Weeks-UNC/shapemapper2/releases/download/2.2.0/shapemapper2-2.2.tar.gz
RUN mkdir -p /opt/shapemapper && \
    curl -fSL "$SM_URL" | tar xz -C /opt/shapemapper --strip-components=1 && \
    find /opt/shapemapper -name '._*' -delete && \
    ln -s /opt/shapemapper/shapemapper /usr/local/bin/shapemapper

WORKDIR /work
ENTRYPOINT []
