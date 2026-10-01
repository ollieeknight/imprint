#!/usr/bin/env nextflow

include { ALIGN }          from './subworkflows/align'
include { QC }             from './subworkflows/qc'
include { SOMATIC_CALLING }  from './subworkflows/somatic_calling'
include { CHARACTERISATION } from './subworkflows/characterisation'
include { ANNOTATION }       from './subworkflows/annotation'
include { PON_GENERATION } from './subworkflows/pon'
include { DUPCALLER } from './subworkflows/dupcaller'
include { EMIT_OUTPUT_MANIFEST; EMIT_PROVENANCE } from './modules/reporting'
include { listParam; effectiveMaxZeroQualFraction; callingBed } from './subworkflows/common'


def preflight_samplesheet(String path) {
    def required = ['patient', 'cell_type', 'status', 'fastq_1', 'fastq_2']
    def parser = new nextflow.util.CsvParser().setSeparator(',').setQuote('"').setStrip(true)
    def lines = new File(path).readLines().findAll { v -> !v.trim().isEmpty() }
    if (lines.isEmpty()) error "Samplesheet is empty: ${path}"

    def headers = parser.parse(lines[0])
    def missing = required - headers
    if (missing) error "Samplesheet missing required columns: ${missing.join(', ')}"

    def rows = lines.tail().withIndex().collect { line, i ->
        def vals = parser.parse(line)
        if (vals.size() != headers.size())
            error "Samplesheet row ${i + 2}: expected ${headers.size()} fields, got ${vals.size()}"
        [headers, vals].transpose().collectEntries()
    }

    if (rows.isEmpty()) error "Samplesheet has headers but no sample rows: ${path}"

    rows.each { row ->
        required.each { col ->
            if (row[col] == null || row[col].trim().isEmpty())
                error "Samplesheet contains an empty '${col}' value"
        }
        ['patient', 'cell_type'].each { col ->
            if (!(row[col] ==~ /[A-Za-z0-9][A-Za-z0-9_-]*/))
                error "Samplesheet ${col} must be letters, digits, '_' or '-', starting with a letter or digit; got '${row[col]}'"
        }
        if (row.status != '0' && row.status != '1')
            error "Samplesheet row for ${row.patient}_${row.cell_type}: status must be '0' (normal) or '1' (tumor), got '${row.status}'"

        ['fastq_1', 'fastq_2'].each { col ->
            if (!file(row[col]).exists())
                error "FASTQ not found for ${row.patient}_${row.cell_type}: ${row[col]}"
        }
    }

    // Also catches repeated rows and fastq_1 == fastq_2.
    def fastq_owners = [:].withDefault { [] }
    rows.each { row ->
        ['fastq_1', 'fastq_2'].each { col ->
            fastq_owners[row[col]] << "${row.patient}_${row.cell_type}:${col}"
        }
    }
    fastq_owners.findAll { _fastqPath, owners -> owners.size() > 1 }.each { fastqPath, owners ->
        error "FASTQ is assigned more than once (${owners.join(', ')}): ${fastqPath}"
    }

    rows.groupBy { r -> "${r.patient}\u0000${r.cell_type}" }.each { key, group ->
        def statuses = group.collect { r -> r.status }.unique()
        if (statuses.size() > 1)
            error "Sample ${key.replace('\u0000', '_')} has conflicting status values: ${statuses.join(', ')}"
    }

    rows.groupBy { r -> r.patient }.each { patient, patRows ->
        def normals = patRows.findAll { r -> r.status == '0' }.collect { r -> r.cell_type }.unique()
        if (normals.size() != 1)
            error "Patient ${patient}: ${normals.size()} normal samples (status=0) found; exactly one required"
    }

    // IDs join names with '_', so different rows can produce the same ID (HC_01 + CD4 vs HC + 01_CD4).
    def samples = rows.collect { r -> [r.patient, r.cell_type, r.status] }.unique()
    def normalOf = samples.findAll { s -> s[2] == '0' }.collectEntries { s -> [s[0], s[1]] }
    def ids = samples.collect { s -> "${s[0]}_${s[1]}" } +
        samples.findAll { s -> s[2] == '1' }.collect { s -> "${s[0]}_${s[1]}_v_${normalOf[s[0]]}" }
    def clashes = ids.groupBy { id -> id.toLowerCase() }.findAll { _id, group -> group.size() > 1 }
    if (clashes)
        error "Samplesheet names produce colliding sample or pair IDs (case-insensitive): ${clashes.values().flatten().unique().join(', ')}"

    rows
}

def asList(value) {
    (value instanceof Collection) ? value : [value]
}

def readFai(String path) {
    file(path).readLines().collectEntries { line ->
        def fields = line.split('\t')
        [(fields[0]): fields[1] as long]
    }
}

def requireFileParam(String name, value) {
    if (value == null || !java.nio.file.Files.isRegularFile(file(value)))
        error "Required file '${name}' not found: ${value}"
}

def requireDirectoryParam(String name, value) {
    if (value == null || !java.nio.file.Files.isDirectory(file(value)))
        error "Required directory '${name}' not found: ${value}"
}

def requireIndexedVcf(String name, value) {
    requireFileParam(name, value)
    if (!file("${value}.tbi").exists())
        error "Index for '${name}' not found: ${value}.tbi"
}

def validateDuplexAllowlist(value) {
    requireFileParam('dupcaller_umi_allowlist', value)
    def umis = file(value).readLines().collect { line -> line.trim().toUpperCase() }.findAll { umi -> umi }
    if (umis.size() != 32 || umis.unique().size() != 32 || umis.any { umi -> !(umi ==~ /[ACGT]{8}/) })
        error "DupCaller allowlist must contain exactly 32 unique 8 bp A/C/G/T UMIs: ${value}"
    def minDistance = umis.withIndex().collectMany { left, i ->
        umis.drop(i + 1).collect { right -> (0..<8).count { k -> left[k] != right[k] } }
    }.min()
    if (minDistance < 3)
        error "DupCaller allowlist minimum Hamming distance must be at least 3; observed ${minDistance}: ${value}"
}

def sampleManifestRecord(meta, String kind, String path) {
    [
        scope: 'sample',
        id: meta.id,
        metadata: [
            sample_id: meta.id,
            donor: meta.donor,
            cell_type: meta.cell_type,
            status: meta.status
        ],
        kind: kind,
        path: path
    ]
}

def pairManifestRecord(meta, String kind, String path) {
    [
        scope: 'pair',
        id: meta.pair_id,
        metadata: [
            pair_id: meta.pair_id,
            donor: meta.donor,
            tumor_id: meta.tumor_id,
            normal_id: meta.normal_id,
            tumor_cell_type: meta.tumor_cell_type,
            normal_cell_type: meta.normal_cell_type
        ],
        kind: kind,
        path: path
    ]
}

def dupcallerManifestRecords(meta, dir, String kindPrefix, String publishRoot) {
    def records = []
    dir.eachFileRecurse(groovy.io.FileType.FILES) { produced ->
        def relative = dir.relativize(produced).toString()
        // SigProfilerPlotting appends the sample name to filenames.
        def suffix = produced.name
            .replaceAll(java.util.regex.Pattern.quote(meta.pair_id), '')
            .toLowerCase()
            .replaceAll('[^a-z0-9]+', '_')
            .replaceAll('^_+|_+$', '')
        records << pairManifestRecord(meta, "${kindPrefix}_${suffix}", "${publishRoot}/${relative}")
    }
    records.sort { record -> record.path }
}

def donorManifestRecord(meta, String kind, String path) {
    [
        scope: 'donor',
        id: meta.donor,
        metadata: [donor: meta.donor],
        kind: kind,
        path: path
    ]
}

// Channel helpers: one manifest record per published file, at its publishDir path.
def sampleRecords(ch, String kind, String subdir) {
    ch.flatMap { meta, files ->
        asList(files).collect { f -> sampleManifestRecord(meta, kind, "${meta.donor}/samples/${meta.cell_type}/${subdir}/${f.name}") }
    }
}

def pairRecords(ch, String kind, String subdir) {
    ch.flatMap { meta, files ->
        asList(files).collect { f -> pairManifestRecord(meta, kind, "${meta.donor}/pairs/${meta.pair_dir}/${subdir}/${f.name}") }
    }
}

// [meta, file, index] recorded as `kind` and `kind_index`.
def indexedPairRecords(ch, String kind, String subdir) {
    pairRecords(ch.map { meta, data, _index -> [meta, data] }, kind, subdir)
        .mix(pairRecords(ch.map { meta, _data, index -> [meta, index] }, "${kind}_index", subdir))
}

def donorRecords(ch, String kind, String subdir) {
    ch.map { meta, f -> donorManifestRecord(meta, kind, "${meta.donor}/${subdir}/${f.name}") }
}

def cohortRecords(ch, String kind, String subdir) {
    ch.flatMap { files ->
        asList(files).collect { f -> [scope: 'cohort', id: params.cohort_name, kind: kind, path: "cohort/${subdir}/${f.name}"] }
    }
}

workflow {

    def resolved_extras = listParam(params.extras).collect { v -> v.toLowerCase() }.unique().sort()
    def unknownExtras = resolved_extras - ['hla', 'kir', 'mixcr', 'pathseq', 'telseq']
    if (unknownExtras)
        error "Unknown --extras value(s): ${unknownExtras.join(', ')}"
    def run_mixcr   = 'mixcr' in resolved_extras
    def run_telseq  = 'telseq' in resolved_extras
    def run_hla     = 'hla' in resolved_extras
    def run_kir     = 'kir' in resolved_extras
    def run_pathseq = 'pathseq' in resolved_extras

    if (params.samplesheet == null) error "Please provide a samplesheet using --samplesheet <path>"
    if (!file(params.samplesheet).exists())  error "Samplesheet not found: ${params.samplesheet}"
    def samplesheet_rows = preflight_samplesheet(params.samplesheet)

    if (params.genome && params.dupcaller)
        error "WGS and DupCaller are not supported together"
    def selectedProfiles = workflow.profile?.tokenize(',')?.collect { profile -> profile.trim() } ?: []
    def selectedKits = selectedProfiles.findAll { profile -> profile in ['agilent_v6', 'agilent_v7', 'agilent_v8', 'twist_v2', 'xgen_exome_v2', 'wgs'] }
    if (!(params.trim_front.toString() ==~ /\d+/))
        error "trim_front must be a non-negative whole number, got: ${params.trim_front}"
    def trimFront = params.trim_front as Integer
    if (trimFront > 0)
        log.info "Trimming ${trimFront} bp from the 5' end of both reads; aligned reads are ${trimFront} bp shorter than sequenced"
    if (params.dupcaller && trimFront != 0)
        error "trim_front applies to the bulk path only; DupCaller mode trims its own UMIs in DUPCALLER_TRIM_LANE"
    if (params.dupcaller && selectedKits != ['xgen_exome_v2'])
        error "DupCaller is restricted to xGen UDSeq; use -profile slurm,xgen_exome_v2,dupcaller"
    def captureParams = ['intervals_bed', 'bait_intervals', 'target_intervals', 'padded_intervals_bed']
        .findAll { name -> params[name] }
    if (params.genome && captureParams)
        error "WGS cannot be combined with capture-kit settings (${captureParams.join(', ')}); use only -profile wgs"

    def effectiveLibraryMode = params.dupcaller ? 'dupcaller' : 'bulk'
    def effectiveMode = params.genome ? 'wgs' : (params.off_target ? 'wes_off_target' : 'wes')
    log.info "imprint: ${effectiveLibraryMode}, ${effectiveMode}${params.kit_name ? ' [' + params.kit_name + ']' : ''}"

    requireFileParam('ref_fasta', params.ref_fasta)
    requireFileParam('genome_fai', params.genome_fai)
    requireFileParam('ref_dict', params.ref_dict)
    if (!file("${params.bwa_mem3_index}.ann").exists())
        error "BWA-MEM3 index files (.ann) not found at prefix: ${params.bwa_mem3_index}"
    if (params.dupcaller) {
        requireFileParam('dupcaller_ref_h5', params.dupcaller_ref_h5)
        requireFileParam('dupcaller_tn_h5', params.dupcaller_tn_h5)
        requireFileParam('dupcaller_hp_h5', params.dupcaller_hp_h5)
        requireFileParam('dupcaller_str_h5', params.dupcaller_str_h5)
        requireFileParam('dupcaller_dbs_h5', params.dupcaller_dbs_h5)
        def expectedH5Names = [
            dupcaller_ref_h5: "${file(params.ref_fasta).name}.ref.h5",
            dupcaller_tn_h5: "${file(params.ref_fasta).name}.tn.h5",
            dupcaller_hp_h5: "${file(params.ref_fasta).name}.hp.h5",
            dupcaller_str_h5: "${file(params.ref_fasta).name}.str.h5",
            dupcaller_dbs_h5: "${file(params.ref_fasta).name}.dbs.h5",
        ]
        expectedH5Names.each { name, expected ->
            if (file(params[name]).name != expected)
                error "${name} must be named ${expected} so DupCaller can resolve it beside the staged FASTA"
        }
        requireIndexedVcf('dupcaller_germline_vcf', params.dupcaller_germline_vcf)
        def dupcallerNoiseMasks = listParam(params.dupcaller_noise_masks)
        dupcallerNoiseMasks.eachWithIndex { mask, index -> requireIndexedVcf("dupcaller_noise_masks[${index}]", mask) }
        if (dupcallerNoiseMasks) log.info "DupCaller noise masks: ${dupcallerNoiseMasks.join(', ')}"
        else log.warn 'DupCaller noise mask is not configured; masked-site filtering is disabled'
        if (params.dupcaller_indel_epon) requireIndexedVcf('dupcaller_indel_epon', params.dupcaller_indel_epon)
        else log.warn 'DupCaller indel ePoN is not configured; enhanced indel panel filtering is disabled'
        validateDuplexAllowlist(params.dupcaller_umi_allowlist)
        if (params.dupcaller_seed != null && !(params.dupcaller_seed.toString() ==~ /\d+/))
            error "dupcaller_seed must be a non-negative whole number, got: ${params.dupcaller_seed}"
    }
    requireIndexedVcf('somalier_sites', params.somalier_sites)
    if (!params.dupcaller) {
        requireIndexedVcf('gnomad_germline_resource_vcf', params.gnomad_germline_resource_vcf)
        requireIndexedVcf('gnomad_pileup_summaries_vcf', params.gnomad_pileup_summaries_vcf)
    }
    requireIndexedVcf('spliceai_snv_vcf', params.spliceai_snv_vcf)
    requireIndexedVcf('gnomad_exomes_vep_vcf', params.gnomad_exomes_vep_vcf)
    requireIndexedVcf('gnomad_genomes_vep_vcf', params.gnomad_genomes_vep_vcf)
    requireIndexedVcf('alphamissense_tsv', params.alphamissense_tsv)
    requireIndexedVcf('dbnsfp_gz', params.dbnsfp_gz)
    requireDirectoryParam('vep_cache', params.vep_cache)
    requireDirectoryParam('vep_plugins_dir', params.vep_plugins_dir)
    requireFileParam('multiqc_config', params.multiqc_config)
    if (!params.dupcaller) {
        requireIndexedVcf('pon_vcf', params.pon_vcf)
        requireIndexedVcf('dbsnp', params.dbsnp)
    }
    if (!params.dupcaller && !params.skip_bqsr) {
        requireIndexedVcf('known_indels_mills', params.known_indels_mills)
        requireIndexedVcf('known_snps_1000g', params.known_snps_1000g)
    }
    if (params.dupcaller && !params.skip_bqsr)
        error "dupcaller mode requires skip_bqsr = true; base recalibration is not permitted before DupCaller"
    if (params.cosmic_vcf) requireIndexedVcf('cosmic_vcf', params.cosmic_vcf)
    if (params.cosmic_noncoding_vcf) requireIndexedVcf('cosmic_noncoding_vcf', params.cosmic_noncoding_vcf)

    def effectiveSpliceAiMode = params.spliceai_indel_vcf ? 'snv_and_indel' : 'snv_only'
    if (params.spliceai_indel_vcf) {
        requireIndexedVcf('spliceai_indel_vcf', params.spliceai_indel_vcf)
        if (file(params.spliceai_indel_vcf).toRealPath() == file(params.spliceai_snv_vcf).toRealPath())
            error "spliceai_indel_vcf points at the SNV VCF; unset it or give the indel VCF"
    } else {
        log.warn "No SpliceAI indel VCF configured; indels get no SpliceAI score"
    }

    if (run_mixcr) requireFileParam('mixcr_license', params.mixcr_license)
    if (run_hla) requireFileParam('hla_reference', params.hla_reference)
    if (run_kir) requireDirectoryParam('kirmapper_db', params.kirmapper_db)

    def ch_pathseq_references = channel.empty()
    if (run_pathseq) {
        def pathseqReferences = [
            pathseq_host_fai    : params.pathseq_host_fai,
            pathseq_host_img    : params.pathseq_host_img,
            pathseq_host_hss    : params.pathseq_host_hss,
            pathseq_microbe_fai : params.pathseq_microbe_fai,
            pathseq_microbe_img : params.pathseq_microbe_img,
            pathseq_microbe_dict: params.pathseq_microbe_dict,
            pathseq_taxonomy    : params.pathseq_taxonomy,
        ]
        pathseqReferences.each { name, referencePath -> requireFileParam(name, referencePath) }

        def primaryContigs = readFai(params.genome_fai)
        def ebvContig = primaryContigs.keySet().find { c -> c == 'chrEBV' || c == 'EBV' }
        if (!ebvContig)
            error "PathSeq requires chrEBV or EBV in the primary alignment reference; neither is present in ${params.genome_fai}"

        def pathseqHostContigs = readFai(params.pathseq_host_fai)
        if (pathseqHostContigs.containsKey(ebvContig))
            error "PathSeq host reference must exclude ${ebvContig}: ${params.pathseq_host_fai}"
        def expectedPathseqHostContigs = primaryContigs.findAll { contig, _length -> contig != ebvContig }
        if (pathseqHostContigs != expectedPathseqHostContigs)
            error "PathSeq host reference must equal the primary alignment reference minus only ${ebvContig}: ${params.pathseq_host_fai}"

        def microbeContigs = readFai(params.pathseq_microbe_fai)
        if (!microbeContigs.containsKey('NC_007605.1'))
            error "PathSeq microbial reference does not contain EBV RefSeq accession NC_007605.1: ${params.pathseq_microbe_fai}"
        if (!microbeContigs.containsKey('NC_006273.2'))
            error "PathSeq microbial reference does not contain CMV RefSeq accession NC_006273.2: ${params.pathseq_microbe_fai}"

        ch_pathseq_references = channel.value([
            file(params.pathseq_host_img),
            file(params.pathseq_host_hss),
            file(params.pathseq_microbe_img),
            file(params.pathseq_microbe_dict),
            file(params.pathseq_taxonomy),
        ])
    }

    if (!params.genome) {
        if (params.intervals_bed == null)
            error "No capture intervals set; add a kit profile from conf/probekits.config, e.g. -profile slurm,xgen_exome_v2"
        def bedParams = ['intervals_bed', 'bait_intervals', 'target_intervals'] + (params.off_target ? ['padded_intervals_bed'] : [])
        bedParams.each { name -> requireFileParam(name, params[name]) }
    }

    if (params.genome && params.verifybamid2_svd_wgs == null)
        error "VerifyBamID2 WGS SVD panel is not configured"

    def donor_sample_counts = samplesheet_rows
        .groupBy { row -> row.patient }
        .collectEntries { donor, rows -> [donor, rows.collect { r -> r.cell_type }.unique().size()] }
    // Preflight enforces one normal per donor.
    def normal_count = donor_sample_counts.size()

    def ch_fastq = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:                 "${row.patient}_${row.cell_type}".toString(),
                donor:              row.patient,
                cell_type:          row.cell_type,
                status:             row.status.toString(),
                donor_sample_count: donor_sample_counts[row.patient],
            ]
            [ meta, file(row.fastq_1), file(row.fastq_2) ]
        }
        .groupTuple(by: 0)
        .flatMap { meta, r1s, r2s ->
            def run_count = r1s.size()
            [ r1s, r2s ].transpose().withIndex().collect { pair, i ->
                def run_meta = meta + [run_count: run_count, run_id: String.format('run%03d', i + 1)]
                [run_meta, pair[0], pair[1]]
            }
        }

    def cohort_dir = file("${params.outdir}/cohort")
    cohort_dir.mkdirs()
    file(params.samplesheet).copyTo(cohort_dir.resolve("samplesheet.csv"))

    def calling_bed = callingBed()

    def effectiveParams = [
        samplesheet             : params.samplesheet,
        outdir                  : params.outdir,
        cohort_name             : params.cohort_name,
        kit_name                : params.kit_name,
        mode                    : effectiveMode,
        library_mode            : effectiveLibraryMode,
        spliceai_mode           : effectiveSpliceAiMode,
        skip_bqsr               : params.dupcaller ? null : params.skip_bqsr,
        calling_intervals       : calling_bed ? calling_bed.toString() : null,
        bait_intervals          : params.genome ? null : params.bait_intervals,
        target_intervals        : params.genome ? null : params.target_intervals,
        scatter_count           : params.dupcaller ? null : (params.genome ? params.genome_scatter_count : params.wes_scatter_count),
        filter_soft_clips       : params.dupcaller ? null : (params.filter_soft_clips == null ? 'autodetect_per_sample' : params.filter_soft_clips),
        optical_dup_dist        : params.optical_dup_dist,
        trim_front              : trimFront,
        extras                  : resolved_extras,
    ].findAll { _name, value -> value != null }
    if (params.dupcaller) {
        effectiveParams += [
            dupcaller_umi_allowlist          : params.dupcaller_umi_allowlist,
            dupcaller_regions                : params.dupcaller_regions,
            dupcaller_max_af                 : params.dupcaller_max_af,
            dupcaller_germline_af_cutoff     : params.dupcaller_germline_af_cutoff,
            dupcaller_min_normal_depth       : params.dupcaller_min_normal_depth,
            dupcaller_max_zero_qual_fraction : params.dupcaller_max_zero_qual_fraction,
            dupcaller_effective_max_zero_qual_fraction: effectiveMaxZeroQualFraction(),
            dupcaller_rescue                 : params.dupcaller_rescue,
            dupcaller_trim_template          : params.dupcaller_trim_template,
            dupcaller_trim_read              : params.dupcaller_trim_read,
            dupcaller_mapq                   : params.dupcaller_mapq,
            dupcaller_window_size            : params.dupcaller_window_size,
            dupcaller_seed                   : params.dupcaller_seed,
        ].findAll { _name, value -> value != null }
    }

    def effectiveReferences = [
        ref_fasta                    : params.ref_fasta,
        genome_fai                  : params.genome_fai,
        ref_dict                    : params.ref_dict,
        bwa_mem3_index              : params.bwa_mem3_index,
        gnomad_germline_resource_vcf: params.dupcaller ? null : params.gnomad_germline_resource_vcf,
        gnomad_pileup_summaries_vcf : params.dupcaller ? null : params.gnomad_pileup_summaries_vcf,
        dbsnp                       : params.dupcaller ? null : params.dbsnp,
        pon_vcf                     : params.dupcaller ? null : params.pon_vcf,
        dupcaller_ref_h5            : params.dupcaller ? params.dupcaller_ref_h5 : null,
        dupcaller_tn_h5             : params.dupcaller ? params.dupcaller_tn_h5 : null,
        dupcaller_hp_h5             : params.dupcaller ? params.dupcaller_hp_h5 : null,
        dupcaller_str_h5            : params.dupcaller ? params.dupcaller_str_h5 : null,
        dupcaller_dbs_h5            : params.dupcaller ? params.dupcaller_dbs_h5 : null,
        dupcaller_germline_vcf      : params.dupcaller ? params.dupcaller_germline_vcf : null,
        dupcaller_noise_masks       : params.dupcaller ? listParam(params.dupcaller_noise_masks) : null,
        dupcaller_indel_epon        : params.dupcaller ? params.dupcaller_indel_epon : null,
        dupcaller_umi_allowlist     : params.dupcaller ? params.dupcaller_umi_allowlist : null,
        known_indels_mills          : params.dupcaller || params.skip_bqsr ? null : params.known_indels_mills,
        known_snps_1000g            : params.dupcaller || params.skip_bqsr ? null : params.known_snps_1000g,
        verifybamid2_svd            : params.genome ? params.verifybamid2_svd_wgs : params.verifybamid2_svd,
        somalier_sites              : params.somalier_sites,
        vep_cache                   : params.vep_cache,
        vep_plugins_dir             : params.vep_plugins_dir,
        alphamissense_tsv           : params.alphamissense_tsv,
        dbnsfp_gz                   : params.dbnsfp_gz,
        spliceai_snv_vcf            : params.spliceai_snv_vcf,
        spliceai_indel_vcf          : params.spliceai_indel_vcf,
        gnomad_exomes_vep_vcf       : params.gnomad_exomes_vep_vcf,
        gnomad_genomes_vep_vcf      : params.gnomad_genomes_vep_vcf,
        cosmic_vcf                  : params.cosmic_vcf,
        cosmic_noncoding_vcf        : params.cosmic_noncoding_vcf,
        hla_reference               : run_hla ? params.hla_reference : null,
        kirmapper_db                : run_kir ? params.kirmapper_db : null,
        pathseq_host_fai            : run_pathseq ? params.pathseq_host_fai : null,
        pathseq_host_img            : run_pathseq ? params.pathseq_host_img : null,
        pathseq_host_hss            : run_pathseq ? params.pathseq_host_hss : null,
        pathseq_microbe_fai         : run_pathseq ? params.pathseq_microbe_fai : null,
        pathseq_microbe_img         : run_pathseq ? params.pathseq_microbe_img : null,
        pathseq_microbe_dict        : run_pathseq ? params.pathseq_microbe_dict : null,
        pathseq_taxonomy            : run_pathseq ? params.pathseq_taxonomy : null,
        mixcr_license               : run_mixcr ? params.mixcr_license : null,
        multiqc_config              : params.multiqc_config,
    ].findAll { _name, value -> value != null }

    def effectiveContainers = [
        align       : params.container_align,
        dupcaller   : params.dupcaller ? params.container_dupcaller : null,
        samtools    : params.container_samtools,
        fastp       : params.container_fastp,
        gatk        : params.container_gatk,
        strelka     : params.dupcaller ? null : params.container_strelka,
        manta       : params.dupcaller ? null : params.container_manta,
        deepsomatic : params.dupcaller ? null : params.container_deepsomatic,
        bcftools    : params.dupcaller ? null : params.container_bcftools,
        vafator     : params.dupcaller ? null : params.container_vafator,
        varlociraptor: params.dupcaller ? null : params.container_varlociraptor,
        vep         : params.container_vep,
        somalier    : params.container_somalier,
        mosdepth    : params.container_mosdepth,
        verifybamid2: params.container_verifybamid2,
        multiqc     : params.container_multiqc,
        yara        : run_hla ? params.container_yara : null,
        optitype    : run_hla ? params.container_optitype : null,
        kir_mapper  : run_kir ? params.container_kir_mapper : null,
        mixcr       : run_mixcr ? params.container_mixcr : null,
        telseq      : run_telseq ? params.container_telseq : null,
    ].findAll { _name, value -> value != null }

    def run_info = [
        pipeline_version: workflow.manifest.version?.toString(),
        commit_id       : workflow.commitId?.toString(),
        repository      : workflow.repository?.toString(),
        revision        : workflow.revision?.toString(),
        profile         : workflow.profile?.toString(),
        command         : workflow.commandLine,
        run_name        : workflow.runName,
        start           : workflow.start.toString(),
        nextflow_version: workflow.nextflow.version.toString(),
        mode            : effectiveMode,
        library_mode    : effectiveLibraryMode,
        extras          : resolved_extras,
        extras_requested: params.extras,
        parameters      : effectiveParams,
        references      : effectiveReferences,
        containers      : effectiveContainers,
    ]
    EMIT_PROVENANCE(run_info)

    ALIGN(ch_fastq)

    def cfg = [
        cohort_name              : params.cohort_name,
        kit_name                 : params.kit_name              ?: 'unknown',
        assay_mode               : effectiveMode,
        calling_intervals        : calling_bed ? calling_bed.name : 'genome-wide',
        off_target               : params.off_target,
        library_mode             : effectiveLibraryMode,
        skip_bqsr                : params.skip_bqsr,
        mosdepth_regions         : calling_bed ? calling_bed.name : 'genome-wide',
        optical_dup_dist         : params.optical_dup_dist,
        trim_front               : trimFront,
    ] + (params.dupcaller ? [
        barcode_chemistry        : 'xgen_udseq_8bp_umi32',
        dupcaller_mapq           : params.dupcaller_mapq,
        dupcaller_min_normal_depth: params.dupcaller_min_normal_depth,
        dupcaller_noise_masks    : listParam(params.dupcaller_noise_masks).collect { mask -> file(mask).name }.join(', ') ?: 'none',
        dupcaller_max_zero_qual_fraction: effectiveMaxZeroQualFraction(),
        dupcaller_rescue         : params.dupcaller_rescue,
    ] : [:])
    def yaml_rows = cfg.collect { k, v -> "    ${k}: '${v}'" }.join('\n')
    def yaml_text = """\
id: 'pipeline_config'
section_name: 'Pipeline Configuration'
description: 'imprint run parameters'
plot_type: 'table'
data:
  Run:
${yaml_rows}
""".stripIndent()

    def ch_config_yaml = channel.of(yaml_text).collectFile(name: 'pipeline_config_mqc.yaml')

    QC(ALIGN.out.analysis_bam, ALIGN.out.reports, ch_config_yaml)

    def ch_sex = QC.out.somalier
        .flatten()
        .filter { f -> f.name.endsWith('.samples.tsv') }
        .splitCsv(header: true, sep: "\t", strip: true)
        .map { row ->
            def clean = row.collectEntries { k, v -> [k.replaceAll('^#', ''), v] }
            def sex
            if (clean.sex == '1') {
                sex = 'XY'
            } else if (clean.sex == '2') {
                sex = 'XX'
            } else {
                def y_depth   = (clean.Y_depth_mean ?: '0').toFloat()
                def mean_depth = (clean.depth_mean   ?: '0').toFloat()
                def y_norm    = (mean_depth > 0) ? y_depth / mean_depth : null
                sex = y_norm == null ? 'UNKNOWN' : (y_norm > 0.1 ? 'XY' : 'XX')
            }
            [ clean.sample_id, sex ]
        }
        .toList()

    // Sorted so the donor merge command, and its -resume hash, is stable.
    def ch_all_bams_per_donor = ALIGN.out.analysis_bam
        .map { meta, bam, bai -> [ groupKey(meta.donor, meta.donor_sample_count as int), meta, bam, bai ] }
        .groupTuple(by: 0)
        .map { donor_key, _metas, bams, bais ->
            def pairs = [bams, bais].transpose().sort { p -> p[0].name }
            [ [donor: donor_key.toString()], pairs.collect { p -> p[0] }, pairs.collect { p -> p[1] } ]
        }

    def branch_bams = ALIGN.out.analysis_bam
        .branch { meta, _bam, _bai ->
            tumor: meta.status == '1'
            normal: meta.status == '0'
        }

    def ch_tumor_enriched = branch_bams.tumor
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
        .join(ALIGN.out.fastp_stats)
        .map { _id, meta, bam, bai, short_inserts, read_length ->
            [meta + [short_inserts: short_inserts, read_length: read_length], bam, bai]
        }

    def ch_paired_bams = ch_tumor_enriched
        .map { meta, bam, bai -> [ meta.donor, meta, bam, bai ] }
        .combine(
            branch_bams.normal.map { meta, bam, bai -> [ meta.donor, meta, bam, bai ] },
            by: 0
        )
        .map { donor, tm, tb, tbai, nm, nb, nbai ->
            def pair_dir = "${tm.cell_type}_v_${nm.cell_type}"
            def paired_meta = [
                pair_id:          "${donor}_${pair_dir}".toString(),
                pair_dir:         pair_dir.toString(),
                donor:            donor,
                tumor_id:         tm.id,
                normal_id:        nm.id,
                tumor_cell_type:  tm.cell_type,
                normal_cell_type: nm.cell_type,
                short_inserts:    tm.short_inserts,
                tumor_count:      tm.donor_sample_count - 1,
            ]
            [ paired_meta, tb, tbai, nb, nbai ]
        }

    CHARACTERISATION(
        ch_all_bams_per_donor,
        ALIGN.out.analysis_bam,
        ALIGN.out.trimmed_reads,
        ALIGN.out.merged_bam,
        ALIGN.out.fastp_stats,
        ch_pathseq_references,
        resolved_extras
    )


    // DupCaller and the bulk callers are mutually exclusive; this is the only place that chooses.
    def ch_mode_records = channel.empty()
    if (params.dupcaller) {
        DUPCALLER(ch_paired_bams)

        ch_mode_records = channel.empty().mix(
            DUPCALLER.out.calls.flatMap { meta, dir -> dupcallerManifestRecords(meta, dir, 'dupcaller_call', "${meta.donor}/pairs/${meta.pair_dir}/dupcaller/calls") },
            DUPCALLER.out.burden.flatMap { meta, dir -> dupcallerManifestRecords(meta, dir, 'dupcaller_burden', "${meta.donor}/pairs/${meta.pair_dir}/dupcaller/burden") },
            DUPCALLER.out.annotated.flatMap { meta, mutation_type, vcf, tbi ->
                def base = "${meta.donor}/pairs/${meta.pair_dir}/dupcaller/annotated"
                [
                    pairManifestRecord(meta, "dupcaller_annotated_${mutation_type}_vcf", "${base}/${vcf.name}"),
                    pairManifestRecord(meta, "dupcaller_annotated_${mutation_type}_vcf_index", "${base}/${tbi.name}")
                ]
            },
            cohortRecords(DUPCALLER.out.cohort_summary, 'dupcaller_summary', 'dupcaller'),
            cohortRecords(DUPCALLER.out.cohort_sbs96, 'dupcaller_sbs96', 'dupcaller')
        )
    } else {
        SOMATIC_CALLING(ch_paired_bams)
        PON_GENERATION(branch_bams.normal, normal_count)
        ANNOTATION(
            SOMATIC_CALLING.out.calling_bams,
            SOMATIC_CALLING.out.mutect2_vcf,
            SOMATIC_CALLING.out.strelka_snv,
            SOMATIC_CALLING.out.strelka_indel,
            SOMATIC_CALLING.out.deepsomatic_vcf,
            ch_sex
        )

        def ch_bqsr_records = SOMATIC_CALLING.out.bqsr_recal.map { meta, role, sample_id, table ->
            def cell_type = role == 'tumor' ? meta.tumor_cell_type : meta.normal_cell_type
            def sample_meta = [id: sample_id, donor: meta.donor, cell_type: cell_type, status: role == 'tumor' ? '1' : '0']
            sampleManifestRecord(sample_meta, 'bqsr_recalibration_table', "${meta.donor}/samples/${cell_type}/qc/bqsr/${table.name}")
        }
        def ch_pon_normal_records = PON_GENERATION.out.normal_vcf.flatMap { meta, vcf, tbi ->
            [
                sampleManifestRecord(meta, 'pon_normal_vcf', "cohort/pon/normals/${vcf.name}"),
                sampleManifestRecord(meta, 'pon_normal_vcf_index', "cohort/pon/normals/${tbi.name}")
            ]
        }

        ch_mode_records = channel.empty().mix(
            ch_bqsr_records,
            ch_pon_normal_records,
            indexedPairRecords(ANNOTATION.out.final_vcf, 'ensemble_vcf', 'variant_calling'),
            indexedPairRecords(ANNOTATION.out.varlociraptor, 'varlociraptor_bcf', 'variant_calling/raw'),
            indexedPairRecords(SOMATIC_CALLING.out.mutect2_raw_vcf, 'mutect2_vcf', 'variant_calling/raw'),
            indexedPairRecords(SOMATIC_CALLING.out.mutect2_vcf, 'mutect2_filtered_vcf', 'variant_calling/raw'),
            pairRecords(SOMATIC_CALLING.out.mutect2_contamination, 'mutect2_contamination_table', 'variant_calling/raw'),
            pairRecords(SOMATIC_CALLING.out.mutect2_segments, 'mutect2_segmentation_table', 'variant_calling/raw'),
            indexedPairRecords(SOMATIC_CALLING.out.strelka_snv.mix(SOMATIC_CALLING.out.strelka_indel), 'strelka2_vcf', 'variant_calling/raw'),
            indexedPairRecords(SOMATIC_CALLING.out.deepsomatic_vcf, 'deepsomatic_vcf', 'variant_calling/raw'),
            indexedPairRecords(SOMATIC_CALLING.out.manta_sv, 'manta_vcf', 'variant_calling'),
            cohortRecords(PON_GENERATION.out.pon.map { vcf, _tbi -> vcf }, 'pon_vcf', 'pon'),
            cohortRecords(PON_GENERATION.out.pon.map { _vcf, tbi -> tbi }, 'pon_vcf_index', 'pon')
        )
    }

    // output_manifest.json records
    def ch_mosdepth_records = QC.out.mosdepth_cov.flatMap { meta, files ->
        asList(files).collect { f ->
            def kind = f.name.endsWith('.mosdepth.summary.txt') ? 'mosdepth_summary' :
                (f.name.endsWith('.mosdepth.region.dist.txt') ? 'mosdepth_region_dist' : 'mosdepth')
            sampleManifestRecord(meta, kind, "${meta.donor}/samples/${meta.cell_type}/qc/mosdepth/${f.name}")
        }
    }

    def ch_riker_records = QC.out.riker_metrics.flatMap { meta, files ->
        asList(files).collect { f ->
            def kind = (f.name.endsWith('.hybcap-metrics.txt') || f.name.endsWith('.wgs-metrics.txt')) ?
                'riker_primary_metrics' : 'riker_metrics'
            sampleManifestRecord(meta, kind, "${meta.donor}/samples/${meta.cell_type}/qc/riker/${f.name}")
        }
    }

    def ch_manifest_records = channel.empty().mix(
        sampleRecords(ALIGN.out.cram.map { meta, cram, _crai -> [meta, cram] }, 'cram', 'alignment'),
        sampleRecords(ALIGN.out.cram.map { meta, _cram, crai -> [meta, crai] }, 'crai', 'alignment'),
        sampleRecords(ALIGN.out.fastp_json, 'fastp_json', 'qc/fastp'),
        sampleRecords(ALIGN.out.fastp_html, 'fastp_html', 'qc/fastp'),
        ch_mosdepth_records,
        sampleRecords(QC.out.selfsm, 'selfsm', 'qc/verifybamid2'),
        ch_riker_records,
        sampleRecords(QC.out.riker_charts, 'riker_chart', 'qc/riker'),
        sampleRecords(QC.out.error_metrics, 'error_rate', 'qc/fgbio'),
        sampleRecords(ALIGN.out.overlap_metrics, 'overlap_metrics', 'qc/fgbio'),
        sampleRecords(ALIGN.out.markdup_metrics,
            params.dupcaller ? 'dupcaller_markdup_metrics' : 'markdup_metrics',
            params.dupcaller ? 'qc/dupcaller' : 'qc/markdup'),
        sampleRecords(ALIGN.out.barcode_metrics, 'dupcaller_barcode_metrics', 'qc/dupcaller'),
        sampleRecords(ALIGN.out.tag_metrics, 'dupcaller_tag_validation', 'qc/dupcaller'),
        sampleRecords(ALIGN.out.cram_tag_metrics, 'dupcaller_cram_tag_validation', 'qc/dupcaller'),
        sampleRecords(CHARACTERISATION.out.pathseq_scores, 'pathseq_scores', 'pathseq'),
        sampleRecords(CHARACTERISATION.out.pathseq_bam, 'pathseq_bam', 'pathseq'),
        sampleRecords(CHARACTERISATION.out.pathseq_filter_metrics, 'pathseq_filter_metrics', 'pathseq'),
        sampleRecords(CHARACTERISATION.out.pathseq_score_metrics, 'pathseq_score_metrics', 'pathseq'),
        sampleRecords(CHARACTERISATION.out.pathseq_score_warnings, 'pathseq_score_warnings', 'pathseq'),
        sampleRecords(CHARACTERISATION.out.telseq, 'telseq', 'telseq'),
        sampleRecords(CHARACTERISATION.out.mixcr_clns, 'mixcr_clns', 'mixcr'),
        sampleRecords(CHARACTERISATION.out.mixcr_clonotypes, 'mixcr_clonotypes', 'mixcr'),
        sampleRecords(CHARACTERISATION.out.mixcr_report, 'mixcr_report', 'mixcr'),
        sampleRecords(CHARACTERISATION.out.mixcr_report_json, 'mixcr_report_json', 'mixcr'),
        sampleRecords(CHARACTERISATION.out.mixcr_step_reports_txt, 'mixcr_step_report', 'mixcr'),
        sampleRecords(CHARACTERISATION.out.mixcr_step_reports_json, 'mixcr_step_report_json', 'mixcr'),
        donorRecords(CHARACTERISATION.out.hla_result, 'optitype_result', 'hla'),
        donorRecords(CHARACTERISATION.out.hla_plot, 'optitype_coverage_plot', 'hla'),
        donorRecords(CHARACTERISATION.out.kir_ncopy, 'kirmapper_ncopy', 'kir'),
        donorRecords(CHARACTERISATION.out.kir_calls, 'kirmapper_calls', 'kir'),
        donorRecords(CHARACTERISATION.out.kir_reports, 'kirmapper_reports', 'kir'),
        donorRecords(CHARACTERISATION.out.kir_raw_archive, 'kirmapper_raw_archive', 'kir'),
        cohortRecords(QC.out.multiqc_report, 'multiqc_report', 'multiqc'),
        cohortRecords(QC.out.multiqc_data, 'multiqc_data', 'multiqc'),
        cohortRecords(QC.out.somalier, 'somalier', 'somalier'),
        ch_mode_records
    ).collect()

    def manifest_info = [
        cohort: params.cohort_name,
        pipeline_version: workflow.manifest.version?.toString(),
        mode: effectiveMode,
        library_mode: effectiveLibraryMode,
        barcode_chemistry: params.dupcaller ? 'xgen_udseq_8bp_umi32' : null
    ]
    EMIT_OUTPUT_MANIFEST(manifest_info, ch_manifest_records)

}
