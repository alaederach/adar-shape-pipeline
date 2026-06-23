// Two-panel skyline control plot.
//   top    — bulk standard-SHAPE (blue) vs sorted/deconvolved SHAPE-ADAR (orange)
//   bottom — ADAR-DMSO editing (red) vs non-ADAR DMSO control (grey), with the
//            sorting cutoff as a horizontal line and every sorted position as a
//            vertical line.
// The DMSO control is the standard-SHAPE sample's Untreated condition (= MOCK-DMSO).
// All profiles are resolved/forwarded by step02; see bin/plot_skyline.py.

process PLOT_SKYLINE {
    tag "${meta.id}"

    publishDir "${params.outdir}/reactivity_profiles/${meta.target}/${meta.cutoff_tag}",
               mode: 'copy'

    input:
    tuple val(meta), path(edited_norm), path(reference_norm), path(adar_dmso), path(dmso_control)

    output:
    tuple val(meta), path("${meta.id}_skyline.png"), emit: plot

    script:
    // candidates mode marks explicit positions (no rate cutoff line); min_rate mode
    // draws the cutoff line and derives the sorted positions from it.
    def cutoff_arg = meta.candidates ? "--sort-positions '${meta.candidates}'"
                   : (meta.min_rate  ? "--cutoff ${meta.min_rate}" : "")
    """
    python3 ${projectDir}/bin/plot_skyline.py \\
        --shape ${reference_norm} \\
        --adar  ${edited_norm} \\
        --adar-dmso ${adar_dmso} \\
        --dmso-control ${dmso_control} \\
        --dmso-control-col Untreated_rate \\
        --refbase ${params.refbase} \\
        ${cutoff_arg} \\
        --title "${meta.target} ${meta.replicate}  (cutoff: ${meta.cutoff_tag})" \\
        --out ${meta.id}_skyline.png
    """
}
