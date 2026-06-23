// Fail-fast validation for params and samplesheet rows, shared by main.nf (step 01)
// and step02.nf. Classes under lib/ are auto-loaded onto the Nextflow classpath.
// These methods do NO file I/O (file existence is checked in the workflow via
// Nextflow's file(..., checkIfExists: true)); they validate columns and values and
// throw IllegalArgumentException with an actionable message.

class Validate {

    static final List<String> NUCS         = ['A', 'C', 'G', 'U', 'T']
    static final List<String> SAMPLE_TYPES  = ['ADAR-SHAPE', 'ADAR-DMSO', 'standard-SHAPE']
    static final List<String> LIBRARY_TYPES = ['targeted-amplicon', 'random-primed']

    private static boolean posInt(Object v) {
        v != null && v.toString().isInteger() && (v.toString() as int) >= 1
    }
    private static boolean nonNegInt(Object v) {
        v != null && v.toString().isInteger() && (v.toString() as int) >= 0
    }
    private static String s(Object v) { v == null ? '' : v.toString().trim() }

    /** Validate pipeline params common to both steps. Call once at workflow start. */
    static void params(Map p) {
        def errs = []
        if (!(s(p.refbase).toUpperCase() in NUCS))
            errs << "refbase must be one of ${NUCS} (got '${p.refbase}')"
        if (!(s(p.readbase).toUpperCase() in NUCS))
            errs << "readbase must be one of ${NUCS} (got '${p.readbase}')"
        if (s(p.refbase) && s(p.refbase).equalsIgnoreCase(s(p.readbase)))
            errs << "refbase and readbase must differ (both '${p.refbase}')"
        if (!posInt(p.min_depth))
            errs << "min_depth must be a positive integer (got '${p.min_depth}')"
        if (!nonNegInt(p.min_mutation_separation))
            errs << "min_mutation_separation must be an integer >= 0 (got '${p.min_mutation_separation}')"
        if (errs)
            throw new IllegalArgumentException(
                "Invalid parameter(s):\n  - " + errs.join("\n  - "))
    }

    /** Validate one step-01 samplesheet row (column presence + values). */
    static void step01Row(Map row) {
        def ctx = rowCtx(row)
        requireColumns(row, ['gene', 'reporter', 'sample_type', 'replicate',
                             'fasta', 'mod_fastq', 'library_type'], ctx)
        def errs = []
        if (!(s(row.sample_type) in SAMPLE_TYPES))
            errs << "sample_type '${row.sample_type}' invalid; expected one of ${SAMPLE_TYPES}"
        if (!(s(row.library_type) in LIBRARY_TYPES))
            errs << "library_type '${row.library_type}' invalid; expected one of ${LIBRARY_TYPES}"
        if (!s(row.mod_fastq))
            errs << "mod_fastq is empty"
        // targeted-amplicon requires primers, else ShapeMapper would silently fall back
        // to random-primer trimming and give wrong results.
        if (s(row.library_type) == 'targeted-amplicon' && !s(row.directed_primers_file))
            errs << "library_type 'targeted-amplicon' requires directed_primers_file " +
                    "(or set library_type to 'random-primed')"
        if (s(row.sample_type) == 'standard-SHAPE' && !s(row.unt_fastq))
            errs << "standard-SHAPE needs an untreated control (unt_fastq is empty)"
        if (errs) throw new IllegalArgumentException("${ctx}:\n  - " + errs.join("\n  - "))
    }

    /** Validate one step-02 samplesheet row (column presence + cutoff values). */
    static void step02Row(Map row) {
        def ctx = rowCtx(row)
        requireColumns(row, ['gene', 'reporter', 'replicate', 'fasta'], ctx)
        def errs = []
        def minRate    = s(row.min_rate)
        def candidates = s(row.candidates)
        if (!minRate && !candidates)
            errs << "provide either min_rate or candidates"
        if (minRate) {
            if (!minRate.isNumber())
                errs << "min_rate must be a number in (0,1] (got '${minRate}')"
            else {
                def r = minRate.toDouble()
                if (r <= 0 || r > 1) errs << "min_rate must be in (0,1] (got ${r})"
            }
        }
        if (candidates && !candidates.split(',').every { it.trim().isInteger() })
            errs << "candidates must be a comma-separated list of integer positions (got '${candidates}')"
        if (errs) throw new IllegalArgumentException("${ctx}:\n  - " + errs.join("\n  - "))
    }

    static void requireColumns(Map row, List<String> required, String ctx) {
        def missing = required.findAll { !row.containsKey(it) || s(row[it]) == '' }
        if (missing)
            throw new IllegalArgumentException(
                "${ctx}: missing/empty required column(s): ${missing.join(', ')}. " +
                "Columns present: ${row.keySet().join(', ')}")
    }

    private static String rowCtx(Map row) {
        "Samplesheet row [${s(row.gene)}_${s(row.reporter)} ${s(row.sample_type)} ${s(row.replicate)}]"
    }
}
