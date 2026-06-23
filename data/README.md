# Demo data — STK4 AluSq2

This folder holds the bundled demo dataset. The reference (`STK4-singleAlus-AluSq2-corrected.fa`)
and primers (`primers/STK4_AluSq2_primers.txt`) are committed; the **FASTQs are not in git**
(they are large binaries). The directory tree under `fastqs/` is kept via `.gitkeep` files so
you can see where each sample belongs.

## Get the FASTQs

```bash
bin/download_data.sh
```

This downloads the dataset from **UNC Dataverse**
([DOI 10.15139/S3/E2BTCM](https://doi.org/10.15139/S3/E2BTCM)) and unpacks the 12 FASTQs
into `fastqs/`.

## ⚠️ This is downsampled data

These FASTQs are **downsampled to ~300,000 read pairs per file** (~178 MB total) so the demo
runs quickly and stays small enough to host. They are representative of the full data (the
deconvolution and normalization reproduce on them), but they are **not** the complete dataset.

The **full-depth sequencing data will be deposited in GEO/SRA when the paper is published**;
this README will be updated with the accession. To rebuild a downsample from full data, see
`bin/make_demo_data.sh`.

## Layout

```
AluSq2/
  STK4-singleAlus-AluSq2-corrected.fa     reference (committed)
  primers/STK4_AluSq2_primers.txt          amplicon primers (committed)
  fastqs/
    ADAR_SHAPE/rep1/   ADAR + SHAPE   (modified, edited reads to sort)
    ADAR_DMSO/rep1/    ADAR + DMSO    (editing-only control)
    MOCK_SHAPE/rep1/   no-ADAR SHAPE  (standard-SHAPE modified)
    MOCK_DMSO/rep1/    no-ADAR DMSO   (standard-SHAPE untreated / DMSO control)
```
