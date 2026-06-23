// ─────────────────────────────────────────────────────────────────────────────
// Rank candidate A->G edit positions from a sample's ShapeMapper profile, sorted
// by editing rate. The published table (especially for the ADAR-DMSO sample) is
// what you inspect to choose the step-02 read-sorting cutoff — a single min_rate
// threshold, or an explicit list of candidate positions.
// ─────────────────────────────────────────────────────────────────────────────

process RANK_POSITIONS {
    tag "${meta.id}"

    publishDir "${params.outdir}/edit_positions", mode: 'copy'

    input:
    tuple val(meta), path(profile)

    output:
    tuple val(meta), path("${meta.id}_ranked_edit_positions.tsv"), emit: ranked

    script:
    """
    python3 ${projectDir}/bin/rank_edit_positions.py \\
        --profile ${profile} \\
        --out     ${meta.id}_ranked_edit_positions.tsv \\
        --refbase ${params.refbase}
    """
}
