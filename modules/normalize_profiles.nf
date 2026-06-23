// Normalize the edited-molecule reactivity profile against the standard-SHAPE
// reference profile (so the edited structure is directly comparable to the bulk
// non-edited structure). tonorm = reference (standard-SHAPE), toscale = edited.
//
// Both outputs come from the SAME normalize call, so they share one normalization
// scale: --normout is the reference re-normalized, --scaleout is the edited scaled
// by that same factor. The skyline control plot consumes both.

process NORMALIZE_PROFILES {
    tag "${meta.id}"

    publishDir "${params.outdir}/reactivity_profiles/${meta.target}/${meta.cutoff_tag}",
               mode: 'copy'

    input:
    tuple val(meta), path(edited_profile), path(reference_profile)

    output:
    tuple val(meta), path("${meta.id}_normalized_profile.txt"), emit: profile
    tuple val(meta),
          path("${meta.id}_normalized_profile.txt"),
          path(ref_out),
          emit: plotdata

    script:
    ref_out = "${meta.target}_${meta.replicate}_${meta.cutoff_tag}_reference_normalized_profile.txt"
    """
    python3 ${projectDir}/bin/normalize_profiles.py \\
        --warn-on-error \\
        --tonorm   ${reference_profile} \\
        --normout  ${ref_out} \\
        --toscale  ${edited_profile} \\
        --scaleout ${meta.id}_normalized_profile.txt
    """
}
