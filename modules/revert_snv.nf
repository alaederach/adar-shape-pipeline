// Revert A->G (readbase->refbase) edits back to A within a region, removing the ADAR
// signature so the re-ShapeMapper step sees only SHAPE adducts, not editing.
// meta.which ('mod'|'unt') keeps the two reverted BAMs distinct.

process REVERT_SNV {
    tag "${meta.id}_${meta.which}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.which}_reverted.bam"),
          path("${meta.id}_${meta.which}_reverted.bam.bai"),
          emit: bam

    script:
    def out = "${meta.id}_${meta.which}_reverted.bam"
    // Region: samplesheet coordinates if given, else the whole reference (from BAM header).
    // The whole-reference fallback is only well-defined for a single-contig reference;
    // if the reference has multiple contigs and no explicit coordinates were given, fail
    // with a clear message rather than silently using the first contig's length.
    def region_resolve = (meta.alu_start != null && meta.alu_end != null)
        ? "RS=${meta.alu_start}; RE=${meta.alu_end}"
        : """NSQ=\$(samtools view -H ${bam} | grep -c '^@SQ' || true)
    if [ "\$NSQ" -ne 1 ]; then
        echo "ERROR: reference has \$NSQ contigs; whole-reference edit reversion is only defined for a single-contig reference. Add alu_start/alu_end to the step-02 samplesheet to set the region explicitly." >&2
        exit 1
    fi
    RS=1
    RE=\$(samtools view -H ${bam} | sed -n 's/.*\\tLN:\\([0-9]*\\).*/\\1/p' | head -1)"""
    """
    ${region_resolve}
    python3 ${projectDir}/bin/pysam_revert_snv_with_qc.py \\
        --input        ${bam} \\
        --out          ${out} \\
        --refbase      ${params.refbase} \\
        --readbase     ${params.readbase} \\
        --region-start \$RS \\
        --region-end   \$RE \\
        --qc-out       ${meta.id}_${meta.which}_revert_qc.tsv
    python3 -c "import pysam, os; pysam.sort('-o','sorted.bam','${out}'); os.rename('sorted.bam','${out}'); pysam.index('${out}')"
    """
}
