#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ═══════════════════════════════════════════════════════════════════════════
// adar-shape-pipeline · STEP 01 — ShapeMapper
//
// Runs ShapeMapper on each sample in the samplesheet and publishes the
// reactivity profiles. For ADAR samples, aligned reads are also emitted so that
// STEP 02 (read sorting) can use them — STEP 02 is run SEPARATELY, after you
// inspect these step-01 profiles and choose an A→G editing-rate cutoff.
// ═══════════════════════════════════════════════════════════════════════════

include { SHAPEMAPPER }    from './modules/shapemapper'
include { RANK_POSITIONS } from './modules/rank_positions'

// Resolve a samplesheet path: keep absolute paths as-is; make relative paths
// relative to the pipeline directory so the pipeline can be launched from anywhere.
def resolvePath = { p -> !p?.trim() ? '' : (p.trim().startsWith('/') ? p.trim() : "${projectDir}/${p.trim()}") }

workflow {

    Validate.params(params)
    if (!params.samplesheet) error "STEP 01 requires --samplesheet <file.csv>"

    Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            Validate.step01Row(row)
            def is_adar = row.sample_type in ['ADAR-SHAPE', 'ADAR-DMSO']
            def meta = [
                id            : "${row.gene}_${row.reporter}_${row.sample_type}_${row.replicate}",
                target        : "${row.gene}_${row.reporter}",
                sample_type   : row.sample_type,
                replicate     : row.replicate,
                library_type  : row.library_type,
                minmutsep     : params.min_mutation_separation,
                // standard-SHAPE is modified+untreated; ADAR samples are modified-only.
                paired_input  : (!is_adar && row.unt_fastq?.trim()) as boolean,
                // ADAR samples emit aligned reads — STEP 02 read sorting consumes them.
                output_aligned: is_adar,
            ]
            // Resolve + check inputs so a bad path fails here with a clear message.
            def fasta   = file(resolvePath(row.fasta), checkIfExists: true)
            def mod_dir = resolvePath(row.mod_fastq)
            if (!file(mod_dir).exists())
                error "${meta.id}: mod_fastq folder not found: ${mod_dir}"
            def unt_dir = resolvePath(row.unt_fastq)
            if (meta.paired_input && !file(unt_dir).exists())
                error "${meta.id}: unt_fastq folder not found: ${unt_dir}"
            def primers = resolvePath(row.directed_primers_file)
            if (primers && !file(primers).exists())
                error "${meta.id}: directed_primers_file not found: ${primers}"
            [ meta, fasta, mod_dir, unt_dir, primers ]
        }
        .set { samples }

    SHAPEMAPPER(samples)

    // Rank A->G edit positions for the ADAR samples. Inspect the published tables
    // (the ADAR-DMSO one especially) to choose the step-02 read-sorting cutoff —
    // either a single min_rate or an explicit list of candidate positions.
    RANK_POSITIONS(
        SHAPEMAPPER.out.profiles
            .filter { meta, profile -> meta.sample_type in ['ADAR-SHAPE', 'ADAR-DMSO'] }
    )
}
