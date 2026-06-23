// Sort reads into "edited" (A->G) within an Alu region. Candidate edit positions are
// defined either by a min_rate threshold on the ADAR-DMSO profile (--min-rate1) or by an
// explicit list of positions (--candidates1) — whichever the step02 samplesheet provides
// (candidates take precedence). One region per construct (single Alu).

process QUERY_SORTED_BAM {
    tag "${meta.id}_${meta.cutoff_tag}"

    publishDir "${params.outdir}/sorted_reads/${meta.target}/${meta.alu_name}/${meta.cutoff_tag}",
               mode: 'copy',
               saveAs: { fn -> (fn.endsWith('.bam') || fn.endsWith('.bam.bai') ||
                                fn.endsWith('_log.txt')) ? fn : null }

    input:
    tuple val(meta),
          path(mod_bam, stageAs: 'mod_sorted.bam'),
          path(mod_bai, stageAs: 'mod_sorted.bam.bai'),
          path(unt_bam, stageAs: 'unt_sorted.bam'),
          path(unt_bai, stageAs: 'unt_sorted.bam.bai'),
          path(mod_profile)

    output:
    tuple val(meta), path("${meta.id}_${meta.cutoff_tag}_edited_mod.bam"),
                     path("${meta.id}_${meta.cutoff_tag}_edited_mod.bam.bai"), emit: edited_mod
    tuple val(meta), path("${meta.id}_${meta.cutoff_tag}_edited_unt.bam"),
                     path("${meta.id}_${meta.cutoff_tag}_edited_unt.bam.bai"), emit: edited_unt
    path("*_log.txt"), optional: true, emit: log

    script:
    // Explicit candidate positions take precedence over the min_rate threshold.
    def cutoff_arg = meta.candidates ? "--candidates1 ${meta.candidates}"
                                     : "--min-rate1 ${meta.min_rate}"
    def no_unedited_arg = params.step02.no_unedited ? "--no-unedited" : ""
    def em = "${meta.id}_${meta.cutoff_tag}_edited_mod.bam"
    def eu = "${meta.id}_${meta.cutoff_tag}_edited_unt.bam"
    // Region: use the samplesheet coordinates if given, else default to the whole
    // reference (length read from the BAM @SQ header).
    def region_resolve = (meta.alu_start != null && meta.alu_end != null)
        ? "START1=${meta.alu_start}; END1=${meta.alu_end}"
        : "START1=1; END1=\$(samtools view -H mod_sorted.bam | sed -n 's/.*\\tLN:\\([0-9]*\\).*/\\1/p' | head -1)"
    """
    ${region_resolve}
    python3 ${projectDir}/bin/query_sorted_bam.py \\
        --modbam     mod_sorted.bam \\
        --untbam     unt_sorted.bam \\
        --modprofile ${mod_profile} \\
        --refbase    ${params.refbase} \\
        --readbase   ${params.readbase} \\
        --edited-min-depth   ${params.step02.edited_min_depth} \\
        --minimum-coverage   ${params.step02.min_coverage_percent} \\
        --max-tolerance      ${params.step02.unedited_max_tolerance} \\
        --unedited-min-depth ${params.step02.unedited_min_depth} \\
        ${no_unedited_arg} \\
        --out-dir . \\
        --start1 \$START1 \\
        --end1   \$END1 \\
        ${cutoff_arg}

    # Edited BAMs are named *_<rate_tag>_atleast..pos_min*kreads.bam (rate_tag is
    # 'minrate*' or 'usercandidates'); rename to predictable output names.
    mv \$(ls mod_sorted_*kreads.bam 2>/dev/null | grep -v _unmodified_ | head -1) ${em}
    mv \$(ls unt_sorted_*kreads.bam 2>/dev/null | grep -v _unmodified_ | head -1) ${eu}
    python3 -c "import pysam; pysam.index('${em}'); pysam.index('${eu}')"
    """
}
