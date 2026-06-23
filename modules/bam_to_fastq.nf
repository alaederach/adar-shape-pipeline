// Convert a coordinate-sorted BAM to a single-end FASTQ (name-sort first, as
// samtools fastq requires). Feeds the reverted reads back into ShapeMapper.

process BAM_TO_FASTQ {
    tag "${meta.id}_${meta.which}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_${meta.which}.fastq"), emit: fastq

    script:
    """
    samtools sort -n -@ ${task.cpus} ${bam} \\
        | samtools fastq -n -@ ${task.cpus} - > ${meta.id}_${meta.which}.fastq
    """
}
