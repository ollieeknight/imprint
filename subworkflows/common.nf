// Functions shared by main.nf and the subworkflows.

def listParam(value) {
    value instanceof Collection
        ? value.findAll { item -> item }
        : (value ? value.toString().tokenize(',').collect { item -> item.trim() }.findAll { item -> item } : [])
}

def effectiveMaxZeroQualFraction() {
    params.dupcaller_max_zero_qual_fraction != null
        ? params.dupcaller_max_zero_qual_fraction
        : (listParam(params.dupcaller_noise_masks) ? 0.5 : 0.1)
}

// References are staged as inputs so that replacing one invalidates the cache.
def fastaRef() {
    [file(params.ref_fasta), file(params.genome_fai), file(params.ref_dict)]
}

def bwaIndex() {
    [file(params.bwa_mem3_index).name, files("${params.bwa_mem3_index}.*")]
}

// [file, index], or [[], []] when the param is unset.
def indexed(path, String suffix = 'tbi') {
    path ? [file(path), file("${path}.${suffix}")] : [[], []]
}

def optionalFile(path) {
    path ? file(path) : []
}

// Padded targets for off-target WES, targets for WES, nothing for WGS.
def callingBed() {
    params.genome ? [] : file(params.off_target ? params.padded_intervals_bed : params.intervals_bed)
}

def tumors(ch_pairs) {
    ch_pairs.map { meta, tb, tbai, _nb, _nbai -> [meta, 'tumor', meta.tumor_id, tb, tbai] }
}

// A normal is shared by every pair in its donor. Emit it once, carrying the first
// pair's meta by name so the task hash is stable.
def normalsOnce(ch_pairs) {
    ch_pairs
        .map { meta, _tb, _tbai, nb, nbai -> [groupKey(meta.normal_id, meta.tumor_count), meta, nb, nbai] }
        .groupTuple(by: 0)
        .map { normal_id, metas, nbs, nbais -> [metas.min { m -> m.pair_id }, 'normal', normal_id.toString(), nbs[0], nbais[0]] }
}

// Rejoin per-sample [meta, role, sample_id, bam, bai] results into [meta, tumor_bam, tumor_bai, normal_bam, normal_bai].
def repair(ch_samples) {
    def split = ch_samples.branch { _meta, role, _sample_id, _bam, _bai ->
        tumor: role == 'tumor'
        normal: role == 'normal'
    }
    split.tumor
        .map { meta, _role, _sample_id, tb, tbai -> [meta.normal_id, meta, tb, tbai] }
        .combine(split.normal.map { _meta, _role, normal_id, nb, nbai -> [normal_id, nb, nbai] }, by: 0)
        .map { _normal_id, meta, tb, tbai, nb, nbai -> [meta, tb, tbai, nb, nbai] }
}
