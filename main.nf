#!/usr/bin/env nextflow

include { ALIGN }          from './subworkflows/align'
include { QC }             from './subworkflows/qc'
include { SOMATIC_CALLING }  from './subworkflows/somatic_calling'
include { CHARACTERISATION } from './subworkflows/characterisation'
include { ANNOTATION }       from './subworkflows/annotation'
include { PON_GENERATION } from './subworkflows/pon'
include { DUPCALLER } from './subworkflows/dupcaller'
include { EMIT_OUTPUT_MANIFEST; EMIT_PROVENANCE } from './modules/reporting'


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
            if (!(row[col] ==~ /[A-Za-z0-9]+/))
                error "Samplesheet ${col} must be strictly alphanumeric ([A-Za-z0-9]+), got '${row[col]}'"
        }
        if (row.status != '0' && row.status != '1')
            error "Samplesheet row for ${row.patient}_${row.cell_type}: status must be '0' (normal) or '1' (tumour), got '${row.status}'"

        if (row.fastq_1 == row.fastq_2)
            error "Samplesheet row for ${row.patient}_${row.cell_type}: fastq_1 and fastq_2 point to the same file"

        ['fastq_1', 'fastq_2'].each { col ->
            if (!file(row[col]).exists())
                error "FASTQ not found for ${row.patient}_${row.cell_type}: ${row[col]}"
        }
    }

    rows.groupBy { r -> "${r.patient}_${r.cell_type}_${r.status}" }.each { key, group ->
        def seen = group.collect { r -> "${r.fastq_1}\u0000${r.fastq_2}" }
        if (seen.size() != seen.unique().size())
            error "Duplicate FASTQ entries for sample ${key}: remove repeated rows"
    }

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
        def normals = patRows
            .findAll { r -> r.status == '0' }
            .collect { r -> "${r.patient}\u0000${r.cell_type}\u0000${r.status}" }
            .unique()
        if (normals.isEmpty())
            error "Patient ${patient}: no normal sample (status=0) found; exactly one required"
        if (normals.size() > 1)
            error "Patient ${patient}: ${normals.size()} normal samples (status=0) found; exactly one required"
    }
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
    if (value == null || !file(value).exists() || !java.nio.file.Files.isRegularFile(file(value)))
        error "Required file '${name}' not found: ${value}"
}

def requireDirectoryParam(String name, value) {
    if (value == null || !file(value).exists() || !java.nio.file.Files.isDirectory(file(value)))
        error "Required directory '${name}' not found: ${value}"
}

def requireIndexedVcf(String name, value) {
    requireFileParam(name, value)
    if (!file("${value}.tbi").exists())
        error "Index for '${name}' not found: ${value}.tbi"
}

def effectiveMaxZeroQualFraction() {
    params.dupcaller_max_zero_qual_fraction != null
        ? params.dupcaller_max_zero_qual_fraction
        : (listParam(params.dupcaller_noise_masks) ? 0.5 : 0.1)
}

def listParam(value) {
    value instanceof Collection
        ? value.findAll { item -> item }
        : (value ? value.toString().tokenize(',').collect { item -> item.trim() }.findAll { item -> item } : [])
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

// A file and its index, recorded as `kind` and `kind_index`.
def indexedPairRecords(meta, String kind, String base, data, index) {
    [
        pairManifestRecord(meta, kind, "${base}/${data.name}"),
        pairManifestRecord(meta, "${kind}_index", "${base}/${index.name}")
    ]
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

workflow {

    def extraNames = [
        mixcr: 'mixcr', telseq: 'telseq', pathseq: 'pathseq',
        hla: 'hla', optitype: 'hla',
        kir: 'kir', kirmapper: 'kir', kir_mapper: 'kir',
    ]
    def requestedExtras = params.extras ? params.extras.toString().tokenize(',').collect { v -> v.trim().toLowerCase() }.findAll { v -> v } : []
    def unknownExtras = requestedExtras.findAll { v -> !extraNames.containsKey(v) }.unique().sort()
    if (unknownExtras)
        error "Unknown --extras value(s): ${unknownExtras.join(', ')}"
    def resolved_extras = requestedExtras.collect { v -> extraNames[v] }.unique().sort()
    def run_mixcr   = 'mixcr' in resolved_extras
    def run_telseq  = 'telseq' in resolved_extras
    def run_hla     = 'hla' in resolved_extras
    def run_kir     = 'kir' in resolved_extras
    def run_pathseq = 'pathseq' in resolved_extras

    def _run_defaults = [
        outdir      : 'imprint_outs',
        extras      : null,
        dupcaller   : false
    ]
    def _overridden = _run_defaults.findAll { k, v -> params[k] != v }
    if (_overridden) {
        def _pad = _overridden.keySet().collect { k -> k.size() }.max()
        log.info "Non-default parameters:\n" + _overridden.collect { k, v ->
            "  --${k.padRight(_pad)}  ${params[k]}  (default: ${v != null ? v : 'null'})"
        }.join("\n")
    }

    if (params.samplesheet == null) error "Please provide a samplesheet using --samplesheet <path>"
    if (!file(params.samplesheet).exists())  error "Samplesheet not found: ${params.samplesheet}"
    preflight_samplesheet(params.samplesheet)

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
    def startupAssayLabel = params.genome ? 'WGS' : (params.off_target ? 'WES off-target' : 'WES target-only')
    log.info "imprint: ${effectiveLibraryMode}, ${startupAssayLabel}${params.kit_name ? ' [' + params.kit_name + ']' : ''}"

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
    requireFileParam('alphamissense_tsv', params.alphamissense_tsv)
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

    def donor_sample_counts = file(params.samplesheet)
        .splitCsv(header: true, strip: true)
        .groupBy { row -> row.patient }
        .collectEntries { donor, rows ->
            [donor, rows.collect { r -> "${r.patient}_${r.cell_type}" }.unique().size()]
        }

    def normal_count = file(params.samplesheet)
        .splitCsv(header: true, strip: true)
        .findAll { row -> row.status == '0' }
        .collect { row -> "${row.patient}_${row.cell_type}" }
        .unique()
        .size()

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

    def effectiveIntervals = params.genome ? null : (params.off_target ? params.padded_intervals_bed : params.intervals_bed)
    def effectiveScatterCount = params.genome ? (params.genome_scatter_count ?: 50) : params.wes_scatter_count

    def effectiveParams = [
        samplesheet             : params.samplesheet,
        outdir                  : params.outdir,
        cohort_name             : params.cohort_name,
        kit_name                : params.kit_name,
        mode                    : effectiveMode,
        library_mode            : effectiveLibraryMode,
        spliceai_mode           : effectiveSpliceAiMode,
        skip_bqsr               : params.dupcaller ? null : params.skip_bqsr,
        calling_intervals       : effectiveIntervals,
        bait_intervals          : params.genome ? null : params.bait_intervals,
        target_intervals        : params.genome ? null : params.target_intervals,
        scatter_count           : params.dupcaller ? null : effectiveScatterCount,
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
        calling_intervals        : effectiveIntervals ? effectiveIntervals.split('/')[-1] : 'genome-wide',
        off_target               : params.off_target,
        library_mode             : effectiveLibraryMode,
        skip_bqsr                : params.skip_bqsr,
        mosdepth_regions         : params.genome ? 'genome-wide' : effectiveIntervals.split('/')[-1],
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

    channel.of(yaml_text)
        .collectFile(name: 'pipeline_config_mqc.yaml')
        .set { ch_config_yaml }

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
            tumour: meta.status == '1'
            normal: meta.status == '0'
        }

    def ch_tumour_enriched = branch_bams.tumour
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
        .join(ALIGN.out.fastp_stats)
        .map { _id, meta, bam, bai, short_inserts, read_length ->
            [meta + [short_inserts: short_inserts, read_length: read_length], bam, bai]
        }

    def ch_paired_bams = ch_tumour_enriched
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
            ]
            [ paired_meta, tb, tbai, nb, nbai ]
        }

    def ch_bulk_paired = params.dupcaller ? channel.empty() : ch_paired_bams
    def ch_dupcaller_paired = params.dupcaller ? ch_paired_bams : channel.empty()
    SOMATIC_CALLING(ch_bulk_paired)
    DUPCALLER(ch_dupcaller_paired, listParam(params.dupcaller_noise_masks), effectiveMaxZeroQualFraction())

    PON_GENERATION(params.dupcaller ? channel.empty() : branch_bams.normal, params.dupcaller ? 0 : normal_count)

    CHARACTERISATION(
        ch_all_bams_per_donor,
        ALIGN.out.analysis_bam,
        ALIGN.out.trimmed_reads,
        ALIGN.out.merged_bam,
        ALIGN.out.fastp_stats,
        ch_pathseq_references,
        resolved_extras
    )

    ANNOTATION(
        // VAFATOR and VarLociraptor use the calling BAMs.
        params.dupcaller ? channel.empty() : SOMATIC_CALLING.out.calling_bams,
        SOMATIC_CALLING.out.mutect2_vcf,
        SOMATIC_CALLING.out.strelka_snv,
        SOMATIC_CALLING.out.strelka_indel,
        SOMATIC_CALLING.out.deepsomatic_vcf,
        ch_sex
    )

    // output_manifest.json records
    def ch_cram_records = ALIGN.out.cram.flatMap { meta, cram, crai ->
        def base = "${meta.donor}/samples/${meta.cell_type}/alignment"
        [
            sampleManifestRecord(meta, 'cram', "${base}/${cram.name}"),
            sampleManifestRecord(meta, 'crai', "${base}/${crai.name}")
        ]
    }

    def ch_fastp_records = ALIGN.out.fastp_json.map { meta, json ->
        sampleManifestRecord(meta, 'fastp_json', "${meta.donor}/samples/${meta.cell_type}/qc/fastp/${json.name}")
    }

    def ch_fastp_html_records = ALIGN.out.fastp_html.map { meta, html ->
        sampleManifestRecord(meta, 'fastp_html', "${meta.donor}/samples/${meta.cell_type}/qc/fastp/${html.name}")
    }

    def ch_mosdepth_records = QC.out.mosdepth_cov.flatMap { meta, files ->
        def base = "${meta.donor}/samples/${meta.cell_type}/qc/mosdepth"
        asList(files).collect { outputFile ->
            def kind = outputFile.name.endsWith('.mosdepth.summary.txt') ? 'mosdepth_summary' :
                (outputFile.name.endsWith('.mosdepth.region.dist.txt') ? 'mosdepth_region_dist' : 'mosdepth')
            sampleManifestRecord(meta, kind, "${base}/${outputFile.name}")
        }
    }

    def ch_selfsm_records = QC.out.selfsm.map { meta, selfsm ->
        sampleManifestRecord(meta, 'selfsm', "${meta.donor}/samples/${meta.cell_type}/qc/verifybamid2/${selfsm.name}")
    }

    def ch_riker_records = QC.out.riker_metrics.flatMap { meta, files ->
        def base = "${meta.donor}/samples/${meta.cell_type}/qc/riker"
        asList(files).collect { outputFile ->
            def kind = (outputFile.name.endsWith('.hybcap-metrics.txt') || outputFile.name.endsWith('.wgs-metrics.txt')) ?
                'riker_primary_metrics' : 'riker_metrics'
            sampleManifestRecord(meta, kind, "${base}/${outputFile.name}")
        }
    }

    def ch_riker_chart_records = QC.out.riker_charts.flatMap { meta, files ->
        asList(files).collect { outputFile ->
            sampleManifestRecord(meta, 'riker_chart', "${meta.donor}/samples/${meta.cell_type}/qc/riker/${outputFile.name}")
        }
    }

    def ch_error_records = QC.out.error_metrics.map { meta, metrics ->
        sampleManifestRecord(meta, 'error_rate', "${meta.donor}/samples/${meta.cell_type}/qc/fgbio/${metrics.name}")
    }

    def ch_overlap_records = ALIGN.out.overlap_metrics.map { meta, metrics ->
        sampleManifestRecord(meta, 'overlap_metrics', "${meta.donor}/samples/${meta.cell_type}/qc/fgbio/${metrics.name}")
    }

    def ch_markdup_records = ALIGN.out.markdup_metrics.map { meta, metrics ->
        def folder = params.dupcaller ? 'dupcaller' : 'markdup'
        def kind = params.dupcaller ? 'dupcaller_markdup_metrics' : 'markdup_metrics'
        sampleManifestRecord(meta, kind, "${meta.donor}/samples/${meta.cell_type}/qc/${folder}/${metrics.name}")
    }

    def ch_dupcaller_barcode_records = ALIGN.out.barcode_metrics.map { meta, metrics ->
        sampleManifestRecord(meta, 'dupcaller_barcode_metrics', "${meta.donor}/samples/${meta.cell_type}/qc/dupcaller/${metrics.name}")
    }

    def ch_dupcaller_tag_records = ALIGN.out.tag_metrics.map { meta, metrics ->
        sampleManifestRecord(meta, 'dupcaller_tag_validation', "${meta.donor}/samples/${meta.cell_type}/qc/dupcaller/${metrics.name}")
    }

    def ch_dupcaller_cram_tag_records = ALIGN.out.cram_tag_metrics.map { meta, metrics ->
        sampleManifestRecord(meta, 'dupcaller_cram_tag_validation', "${meta.donor}/samples/${meta.cell_type}/qc/dupcaller/${metrics.name}")
    }

    def ch_dupcaller_call_records = DUPCALLER.out.calls.flatMap { meta, dir ->
        dupcallerManifestRecords(meta, dir, 'dupcaller_call',
            "${meta.donor}/pairs/${meta.pair_dir}/dupcaller/calls")
    }

    def ch_dupcaller_burden_records = DUPCALLER.out.burden.flatMap { meta, dir ->
        dupcallerManifestRecords(meta, dir, 'dupcaller_burden',
            "${meta.donor}/pairs/${meta.pair_dir}/dupcaller/burden")
    }

    def ch_dupcaller_annotated_records = DUPCALLER.out.annotated.flatMap { meta, mutation_type, vcf, tbi ->
        indexedPairRecords(meta, "dupcaller_annotated_${mutation_type}_vcf", "${meta.donor}/pairs/${meta.pair_dir}/dupcaller/annotated", vcf, tbi)
    }

    def ch_bqsr_records = SOMATIC_CALLING.out.bqsr_recal.map { meta, role, sample_id, table ->
        def cell_type = role == 'tumour' ? meta.tumor_cell_type : meta.normal_cell_type
        def sample_meta = [
            id: sample_id, donor: meta.donor, cell_type: cell_type,
            status: role == 'tumour' ? '1' : '0'
        ]
        sampleManifestRecord(sample_meta, 'bqsr_recalibration_table', "${meta.donor}/samples/${cell_type}/qc/bqsr/${table.name}")
    }

    def ch_pathseq_records = CHARACTERISATION.out.pathseq_scores.map { meta, scores ->
        sampleManifestRecord(meta, 'pathseq_scores', "${meta.donor}/samples/${meta.cell_type}/pathseq/${scores.name}")
    }

    def ch_pathseq_bam_records = CHARACTERISATION.out.pathseq_bam.map { meta, bam ->
        sampleManifestRecord(meta, 'pathseq_bam', "${meta.donor}/samples/${meta.cell_type}/pathseq/${bam.name}")
    }

    def ch_pathseq_filter_records = CHARACTERISATION.out.pathseq_filter_metrics.map { meta, metrics ->
        sampleManifestRecord(meta, 'pathseq_filter_metrics', "${meta.donor}/samples/${meta.cell_type}/pathseq/${metrics.name}")
    }

    def ch_pathseq_score_records = CHARACTERISATION.out.pathseq_score_metrics.map { meta, metrics ->
        sampleManifestRecord(meta, 'pathseq_score_metrics', "${meta.donor}/samples/${meta.cell_type}/pathseq/${metrics.name}")
    }

    def ch_pathseq_warning_records = CHARACTERISATION.out.pathseq_score_warnings.map { meta, warnings ->
        sampleManifestRecord(meta, 'pathseq_score_warnings', "${meta.donor}/samples/${meta.cell_type}/pathseq/${warnings.name}")
    }

    def ch_telseq_records = CHARACTERISATION.out.telseq.map { meta, result ->
        sampleManifestRecord(meta, 'telseq', "${meta.donor}/samples/${meta.cell_type}/telseq/${result.name}")
    }

    def ch_mixcr_records = CHARACTERISATION.out.mixcr_clns.map { meta, clns ->
        sampleManifestRecord(meta, 'mixcr_clns', "${meta.donor}/samples/${meta.cell_type}/mixcr/${clns.name}")
    }.mix(CHARACTERISATION.out.mixcr_clonotypes.flatMap { meta, files ->
            asList(files).collect { outputFile ->
                sampleManifestRecord(meta, 'mixcr_clonotypes', "${meta.donor}/samples/${meta.cell_type}/mixcr/${outputFile.name}")
            }
        })
        .mix(CHARACTERISATION.out.mixcr_report.map { meta, result ->
            sampleManifestRecord(meta, 'mixcr_report', "${meta.donor}/samples/${meta.cell_type}/mixcr/${result.name}")
        })
        .mix(CHARACTERISATION.out.mixcr_report_json.map { meta, result ->
            sampleManifestRecord(meta, 'mixcr_report_json', "${meta.donor}/samples/${meta.cell_type}/mixcr/${result.name}")
        })
        .mix(CHARACTERISATION.out.mixcr_step_reports_txt.flatMap { meta, files ->
            asList(files).collect { outputFile ->
                sampleManifestRecord(meta, 'mixcr_step_report', "${meta.donor}/samples/${meta.cell_type}/mixcr/${outputFile.name}")
            }
        })
        .mix(CHARACTERISATION.out.mixcr_step_reports_json.flatMap { meta, files ->
            asList(files).collect { outputFile ->
                sampleManifestRecord(meta, 'mixcr_step_report_json', "${meta.donor}/samples/${meta.cell_type}/mixcr/${outputFile.name}")
            }
        })

    def ch_hla_records = CHARACTERISATION.out.hla_result.map { meta, result ->
        donorManifestRecord(meta, 'optitype_result', "${meta.donor}/hla/${result.name}")
    }

    def ch_hla_plot_records = CHARACTERISATION.out.hla_plot.map { meta, plot ->
        donorManifestRecord(meta, 'optitype_coverage_plot', "${meta.donor}/hla/${plot.name}")
    }

    def ch_kir_records = CHARACTERISATION.out.kir_ncopy.map { meta, outputFile ->
        donorManifestRecord(meta, 'kirmapper_ncopy', "${meta.donor}/kir/${outputFile.name}")
    }.mix(CHARACTERISATION.out.kir_calls.map { meta, outputFile ->
        donorManifestRecord(meta, 'kirmapper_calls', "${meta.donor}/kir/${outputFile.name}")
    }).mix(CHARACTERISATION.out.kir_reports.map { meta, outputFile ->
        donorManifestRecord(meta, 'kirmapper_reports', "${meta.donor}/kir/${outputFile.name}")
    }).mix(CHARACTERISATION.out.kir_raw_archive.map { meta, outputFile ->
        donorManifestRecord(meta, 'kirmapper_raw_archive', "${meta.donor}/kir/${outputFile.name}")
    })

    def ch_final_call_records = ANNOTATION.out.final_vcf.flatMap { meta, vcf, tbi ->
        indexedPairRecords(meta, 'ensemble_vcf', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling", vcf, tbi)
    }

    def ch_varlociraptor_records = ANNOTATION.out.varlociraptor.flatMap { meta, bcf, csi ->
        indexedPairRecords(meta, 'varlociraptor_bcf', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw", bcf, csi)
    }

    def ch_mutect2_records = SOMATIC_CALLING.out.mutect2_raw_vcf.flatMap { meta, vcf, tbi ->
        indexedPairRecords(meta, 'mutect2_vcf', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw", vcf, tbi)
    }

    def ch_mutect2_filtered_records = SOMATIC_CALLING.out.mutect2_vcf.flatMap { meta, vcf, tbi ->
        indexedPairRecords(meta, 'mutect2_filtered_vcf', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw", vcf, tbi)
    }

    def ch_mutect2_contamination_records = SOMATIC_CALLING.out.mutect2_contamination.map { meta, table ->
        pairManifestRecord(meta, 'mutect2_contamination_table', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw/${table.name}")
    }.mix(SOMATIC_CALLING.out.mutect2_segments.map { meta, table ->
        pairManifestRecord(meta, 'mutect2_segmentation_table', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw/${table.name}")
    })

    def ch_pon_normal_records = PON_GENERATION.out.normal_vcf.flatMap { meta, vcf, tbi ->
        def base = 'cohort/pon/normals'
        [
            sampleManifestRecord(meta, 'pon_normal_vcf', "${base}/${vcf.name}"),
            sampleManifestRecord(meta, 'pon_normal_vcf_index', "${base}/${tbi.name}")
        ]
    }

    def ch_strelka_records = SOMATIC_CALLING.out.strelka_snv
        .mix(SOMATIC_CALLING.out.strelka_indel)
        .flatMap { meta, vcf, tbi ->
            indexedPairRecords(meta, 'strelka2_vcf', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw", vcf, tbi)
        }

    def ch_deepsomatic_records = SOMATIC_CALLING.out.deepsomatic_vcf.flatMap { meta, vcf, tbi ->
        indexedPairRecords(meta, 'deepsomatic_vcf', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling/raw", vcf, tbi)
    }

    def ch_manta_records = SOMATIC_CALLING.out.manta_sv.flatMap { meta, vcf, tbi ->
        indexedPairRecords(meta, 'manta_vcf', "${meta.donor}/pairs/${meta.pair_dir}/variant_calling", vcf, tbi)
    }

    def ch_cohort_records = QC.out.multiqc_report
        .map { report -> [scope: 'cohort', id: params.cohort_name, kind: 'multiqc_report', path: "cohort/multiqc/${report.name}"] }
        .mix(QC.out.multiqc_data.map { data -> [scope: 'cohort', id: params.cohort_name, kind: 'multiqc_data', path: "cohort/multiqc/${data.name}"] })
        .mix(QC.out.somalier.flatMap { files ->
            asList(files).collect { outputFile ->
                [scope: 'cohort', id: params.cohort_name, kind: 'somalier', path: "cohort/somalier/${outputFile.name}"]
            }
        })
        .mix(PON_GENERATION.out.pon.flatMap { _cohort, vcf, tbi ->
            [
                [scope: 'cohort', id: params.cohort_name, kind: 'pon_vcf', path: "cohort/pon/${vcf.name}"],
                [scope: 'cohort', id: params.cohort_name, kind: 'pon_vcf_index', path: "cohort/pon/${tbi.name}"]
            ]
        })
        .mix(DUPCALLER.out.cohort_summary.map { outputFile ->
            [scope: 'cohort', id: params.cohort_name, kind: 'dupcaller_summary', path: "cohort/dupcaller/${outputFile.name}"]
        })
        .mix(DUPCALLER.out.cohort_sbs96.flatMap { files ->
            asList(files).collect { outputFile ->
                [scope: 'cohort', id: params.cohort_name, kind: 'dupcaller_sbs96', path: "cohort/dupcaller/${outputFile.name}"]
            }
        })

    def ch_manifest_records = ch_cram_records
        .mix(ch_fastp_records)
        .mix(ch_fastp_html_records)
        .mix(ch_mosdepth_records)
        .mix(ch_selfsm_records)
        .mix(ch_riker_records)
        .mix(ch_riker_chart_records)
        .mix(ch_error_records)
        .mix(ch_overlap_records)
        .mix(ch_markdup_records)
        .mix(ch_dupcaller_barcode_records)
        .mix(ch_dupcaller_tag_records)
        .mix(ch_dupcaller_cram_tag_records)
        .mix(ch_dupcaller_call_records)
        .mix(ch_dupcaller_burden_records)
        .mix(ch_dupcaller_annotated_records)
        .mix(ch_bqsr_records)
        .mix(ch_pathseq_records)
        .mix(ch_pathseq_bam_records)
        .mix(ch_pathseq_filter_records)
        .mix(ch_pathseq_score_records)
        .mix(ch_pathseq_warning_records)
        .mix(ch_telseq_records)
        .mix(ch_mixcr_records)
        .mix(ch_hla_records)
        .mix(ch_hla_plot_records)
        .mix(ch_kir_records)
        .mix(ch_final_call_records)
        .mix(ch_varlociraptor_records)
        .mix(ch_mutect2_records)
        .mix(ch_mutect2_filtered_records)
        .mix(ch_mutect2_contamination_records)
        .mix(ch_pon_normal_records)
        .mix(ch_strelka_records)
        .mix(ch_deepsomatic_records)
        .mix(ch_manta_records)
        .mix(ch_cohort_records)
        .collect()

    def manifest_info = [
        cohort: params.cohort_name,
        pipeline_version: workflow.manifest.version?.toString(),
        mode: effectiveMode,
        library_mode: effectiveLibraryMode,
        barcode_chemistry: params.dupcaller ? 'xgen_udseq_8bp_umi32' : null
    ]
    EMIT_OUTPUT_MANIFEST(manifest_info, ch_manifest_records)

}
