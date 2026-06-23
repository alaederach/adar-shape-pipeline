// Filter a ShapeMapper aligned SAM to uniquely-mapped reads -> coordinate-sorted,
// indexed BAM. Step02 runs this on both the ADAR-SHAPE (mod) and ADAR-DMSO (unt)
// aligned reads; meta.which ('mod'|'unt') keeps the two outputs distinct.

process SAM_TO_UNIQUE_BAM {
    tag "${meta.id}_${meta.which}"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.which}.unique.sorted.bam"),
          path("${meta.id}_${meta.which}.unique.sorted.bam.bai"),
          emit: bam

    script:
    """
    python3 ${projectDir}/bin/sam_to_unique_sorted_bam.py \\
        --sam     ${sam} \\
        --out     ${meta.id}_${meta.which}.unique.sorted.bam \\
        --threads ${task.cpus}
    """
}
