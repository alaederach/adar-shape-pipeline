// ─────────────────────────────────────────────────────────────────────────────
// ShapeMapper — run on a single sample; publish its reactivity profile(s).
//
// meta keys: id, target, sample_type, replicate, library_type, minmutsep,
//            paired_input (bool), output_aligned (bool)
//
// ShapeMapper 2.2 is x86-64 Linux only. `params.shapemapper_bin` points at the
// executable: a native install on Linux, or a container wrapper elsewhere
// (see the `docker` profile and the README).
// ─────────────────────────────────────────────────────────────────────────────

process SHAPEMAPPER {
    tag "${meta.id}"

    publishDir "${params.outdir}/shapemapper/${meta.id}",
               mode: 'copy',
               saveAs: { fn ->
                   if (fn.endsWith('_profile.txt') || fn.endsWith('.shape') ||
                       fn.endsWith('_mutation_counts.txt')) return fn
                   // ADAR samples: also publish the aligned reads — STEP 02 read sorting
                   // resolves these by convention ({outdir}/shapemapper/{id}/{id}_aligned.sam).
                   if (meta.output_aligned && fn == "${meta.id}_aligned.sam") return fn
                   return null
               }

    input:
    tuple val(meta), path(fasta), val(mod_fastq), val(unt_fastq), val(primers)

    output:
    tuple val(meta), path("*_profile.txt"),               emit: profiles
    tuple val(meta), path("*_mutation_counts.txt"),       emit: counts
    tuple val(meta), path("${meta.id}_aligned.sam"),      emit: aligned
    tuple val(meta), path("*"),                           emit: all

    script:
    // fastq_mode 'folder' (default, step01 raw fastqs) or 'single_end' (step02
    // reverted reads as a single FASTQ via --U).
    def src_flag   = (meta.fastq_mode == 'single_end') ? '--U' : '--folder'
    def mod_arg    = "--modified ${src_flag} ${mod_fastq}"
    def unt_arg    = meta.paired_input ? "--untreated ${src_flag} ${unt_fastq}" : ''
    def primer_arg = (meta.library_type == 'targeted-amplicon' && primers)
        ? "--primers ${primers} --amplicon"
        : "--random-primer-len ${params.random_primer_len}"
    def aligned_arg = meta.output_aligned ? '--output-aligned-reads' : ''

    // ShapeMapper 2.2 names the --output-aligned-reads SAM `{name}_Modified_aligned.sam`,
    // not `{name}_aligned.sam`. Rename it to the declared output. When aligned reads
    // are not requested, the tuple output still requires the file, so touch a placeholder.
    def finalize_sam = meta.output_aligned
        ? "mv \$(ls *_Modified_aligned.sam | head -1) ${meta.id}_aligned.sam"
        : "touch ${meta.id}_aligned.sam"

    """
    ${params.shapemapper_bin} \\
        --name ${meta.id} \\
        --target ${fasta} \\
        ${mod_arg} \\
        ${unt_arg} \\
        --min-mutation-separation ${meta.minmutsep} \\
        ${primer_arg} \\
        --min-depth ${params.min_depth} \\
        --nproc ${task.cpus} \\
        --output-counted-mutations --overwrite \\
        ${aligned_arg} \\
        --out .

    rm -rf shapemapper_temp
    ${finalize_sam}
    """
}
