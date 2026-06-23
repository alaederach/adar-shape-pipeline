#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ═══════════════════════════════════════════════════════════════════════════
// adar-shape-pipeline · STEP 02 — sort + deconvolve → edited-molecule reactivity
//
// Run AFTER step01, once you've inspected the ADAR-DMSO ranked edit table and
// chosen a cutoff. Resolves step01 outputs by convention from --step01_outdir,
// then for each construct:
//   a) SAM -> unique BAM (ADAR-SHAPE mod + ADAR-DMSO unt)
//   b) sort reads into the edited (A->G) set using the chosen cutoff
//   c) revert the A->G edits, re-run ShapeMapper (mod+unt) on the reverted reads
//   d) normalize the resulting reactivity profile against the standard-SHAPE reference
//
//   nextflow run step02.nf -profile <profile> \
//       --samplesheet input/samplesheet_step02.csv --step01_outdir results
// ═══════════════════════════════════════════════════════════════════════════

include { SAM_TO_UNIQUE_BAM as SAM_TO_BAM_MOD } from './modules/sam_to_bam'
include { SAM_TO_UNIQUE_BAM as SAM_TO_BAM_UNT } from './modules/sam_to_bam'
include { QUERY_SORTED_BAM }                    from './modules/query_sorted_bam'
include { REVERT_SNV as REVERT_MOD }            from './modules/revert_snv'
include { REVERT_SNV as REVERT_UNT }            from './modules/revert_snv'
include { BAM_TO_FASTQ as BAM2FQ_MOD }          from './modules/bam_to_fastq'
include { BAM_TO_FASTQ as BAM2FQ_UNT }          from './modules/bam_to_fastq'
include { SHAPEMAPPER as RESHAPE }              from './modules/shapemapper'
include { NORMALIZE_PROFILES }                  from './modules/normalize_profiles'
include { PLOT_SKYLINE }                        from './modules/plot_skyline'

def resolvePath = { p -> !p?.trim() ? null : (p.trim().startsWith('/') ? p.trim() : "${projectDir}/${p.trim()}") }

workflow {

    Validate.params(params)
    if (!params.samplesheet) error "STEP 02 requires --samplesheet <file.csv>"
    if (!params.step01_outdir)
        error "STEP 02 requires --step01_outdir (path to step01's results/ directory)."
    def s01 = resolvePath(params.step01_outdir)

    Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true, quote: '"')
        .map { row ->
            Validate.step02Row(row)
            def base    = "${row.gene}_${row.reporter}"
            def repl    = row.replicate
            def shapeId = "${base}_ADAR-SHAPE_${repl}"
            def dmsoId  = "${base}_ADAR-DMSO_${repl}"
            def stdId   = "${base}_standard-SHAPE_${repl}"

            // Resolve step01 outputs by convention.
            def modSam = file("${s01}/shapemapper/${shapeId}/${shapeId}_aligned.sam", checkIfExists: true)
            def untSam = file("${s01}/shapemapper/${dmsoId}/${dmsoId}_aligned.sam",  checkIfExists: true)
            def dmsoP  = files("${s01}/shapemapper/${dmsoId}/*_profile.txt")   // candidate selection
            def stdP   = files("${s01}/shapemapper/${stdId}/*_profile.txt")    // normalization reference
            if (!dmsoP) error "No ADAR-DMSO profile for ${dmsoId} under ${s01}"
            if (!stdP)  error "No standard-SHAPE profile for ${stdId} under ${s01}"
            def fasta  = file(resolvePath(row.fasta), checkIfExists: true)

            // Cutoff: candidates (explicit positions) take precedence over min_rate.
            def candidates = row.candidates?.trim() ?: null
            def min_rate   = row.min_rate?.trim()   ?: null
            if (!candidates && !min_rate)
                error "Row ${base} ${repl}: provide min_rate or candidates."
            def cutoff_tag = candidates ? "candidates"
                : "minrate${((min_rate.toFloat() * 100) as int).toString().padLeft(2, '0')}"

            // Region optional; whole reference if omitted.
            def alu_name  = row.alu_name?.trim()  ?: 'full'
            def alu_start = row.alu_start?.trim() ? row.alu_start.toInteger() : null
            def alu_end   = row.alu_end?.trim()   ? row.alu_end.toInteger()   : null

            def meta = [ id: "${base}_${repl}_${alu_name}", target: base, replicate: repl,
                         alu_name: alu_name, alu_start: alu_start, alu_end: alu_end,
                         min_rate: min_rate, candidates: candidates, cutoff_tag: cutoff_tag ]
            [ meta, modSam, untSam, dmsoP.first(), stdP.first(), fasta ]
        }
        .set { samples }

    // construct+cutoff key, to re-join the fasta and the reference profile downstream
    def ckey = { m -> "${m.target}_${m.replicate}_${m.cutoff_tag}" }

    // ── (a) SAM -> unique BAM ────────────────────────────────────────────────
    SAM_TO_BAM_MOD( samples.map { m, mod, unt, dp, sp, fa -> [ m + [which: 'mod'], mod ] } )
    SAM_TO_BAM_UNT( samples.map { m, mod, unt, dp, sp, fa -> [ m + [which: 'unt'], unt ] } )

    // ── (b) sort reads into the edited set ───────────────────────────────────
    QUERY_SORTED_BAM(
        SAM_TO_BAM_MOD.out.bam.map { m, b, i -> [ m.id, m, b, i ] }
            .join( SAM_TO_BAM_UNT.out.bam.map { m, b, i -> [ m.id, b, i ] } )
            .join( samples.map { m, mod, unt, dp, sp, fa -> [ m.id, dp ] } )
            .map { id, m, mb, mi, ub, ui, dp -> [ m.subMap(m.keySet() - ['which']), mb, mi, ub, ui, dp ] }
    )

    // ── (c) revert edits → FASTQ → re-ShapeMapper (mod+unt, single-end) ───────
    REVERT_MOD( QUERY_SORTED_BAM.out.edited_mod.map { m, b, i -> [ m + [which: 'mod'], b, i ] } )
    REVERT_UNT( QUERY_SORTED_BAM.out.edited_unt.map { m, b, i -> [ m + [which: 'unt'], b, i ] } )
    BAM2FQ_MOD( REVERT_MOD.out.bam )
    BAM2FQ_UNT( REVERT_UNT.out.bam )

    aux = samples.map { m, mod, unt, dp, sp, fa -> [ ckey(m), fa, sp ] }

    RESHAPE(
        BAM2FQ_MOD.out.fastq.map { m, fq -> [ ckey(m), m, fq ] }
            .join( BAM2FQ_UNT.out.fastq.map { m, fq -> [ ckey(m), fq ] } )
            .join( aux )
            .map { key, m, modfq, untfq, fa, sp ->
                def rm = [ id: "${m.target}_${m.replicate}_${m.cutoff_tag}_edited",
                           target: m.target, replicate: m.replicate, cutoff_tag: m.cutoff_tag,
                           min_rate: m.min_rate, candidates: m.candidates,
                           library_type: 'random-primed', fastq_mode: 'single_end',
                           paired_input: true, output_aligned: false,
                           minmutsep: params.min_mutation_separation ]
                [ rm, fa, modfq.toString(), untfq.toString(), '' ]
            }
    )

    // ── (d) normalize edited reactivity vs standard-SHAPE reference ───────────
    NORMALIZE_PROFILES(
        RESHAPE.out.profiles.map { m, prof -> [ ckey(m), m, prof ] }
            .join( samples.map { m, mod, unt, dp, sp, fa -> [ ckey(m), sp ] } )
            .map { key, m, prof, sp -> [ m, prof, sp ] }
    )

    // ── (e) two-panel skyline control plot ───────────────────────────────────
    //   top    : bulk SHAPE (blue) vs edited SHAPE-ADAR (orange)
    //   bottom : ADAR-DMSO editing (red) vs DMSO control (grey) + cutoff/sort lines.
    // dp = ADAR-DMSO profile (editing); sp = standard-SHAPE profile (its Untreated
    // condition is the non-ADAR DMSO control). Joined back from `samples` by ckey.
    PLOT_SKYLINE(
        NORMALIZE_PROFILES.out.plotdata.map { m, ed, ref -> [ ckey(m), m, ed, ref ] }
            .join( samples.map { m, mod, unt, dp, sp, fa -> [ ckey(m), dp, sp ] } )
            .map { key, m, ed, ref, dp, sp -> [ m, ed, ref, dp, sp ] }
    )
}
