def jsonSafe(value) {
    if (value == null || value instanceof Number || value instanceof Boolean || value instanceof CharSequence) return value
    if (value instanceof Map) {
        return value.collectEntries { key, item ->
            def name = key.toString()
            def redacted = name ==~ /(?i).*(password|passwd|token|secret|api[_-]?key|access[_-]?key).*/
            [(name): redacted && item != null ? '[REDACTED]' : jsonSafe(item)]
        }
    }
    if (value instanceof Collection) return value.collect { item -> jsonSafe(item) }
    if (value.getClass().isArray()) return value.toList().collect { item -> jsonSafe(item) }
    // Nextflow params commonly contain Path, Duration, MemoryUnit and other
    // configuration value types that JsonOutput cannot serialise directly.
    value.toString()
}

def groupManifestRecords(List records, String scope) {
    def grouped = records
        .findAll { record -> record.scope == scope }
        .groupBy { record -> record.id }
        .collect { id, scopedRecords ->
            def metadataValues = scopedRecords.collect { record -> record.metadata ?: [:] }
            def metadataKeys = metadataValues.collectMany { item -> item.keySet() }.unique()
            def metadata = metadataKeys.sort().collectEntries { key ->
                def values = metadataValues.collect { item -> item[key] }.findAll { value -> value != null }.unique()
                if (values.size() > 1) {
                    throw new IllegalArgumentException("Conflicting ${scope} metadata for '${id}', field '${key}': ${values}")
                }
                [(key): values ? values[0] : null]
            }
            def outputs = new TreeMap(scopedRecords
                .groupBy { record -> record.kind }
                .collectEntries { kind, outputRecords ->
                    [kind, outputRecords.collect { record -> record.path }.findAll { path -> path }.unique().sort()]
                })
            metadata + [outputs: outputs]
        }
    grouped.sort { left, right ->
        (left.sample_id ?: left.pair_id ?: left.donor ?: '') <=> (right.sample_id ?: right.pair_id ?: right.donor ?: '')
    }
}

process EMIT_OUTPUT_MANIFEST {
    label 'process_single'
    tag 'cohort'
    publishDir "${params.outdir}/cohort", mode: 'copy'
    executor 'local'

    input:
    val manifest_info
    val records

    output:
    path 'output_manifest.json', emit: output_manifest

    exec:
    def recordList = (records ?: []).flatten().findAll { record -> record }.sort { left, right ->
        "${left.scope ?: ''}\u0000${left.id ?: ''}\u0000${left.kind ?: ''}\u0000${left.path ?: ''}" <=>
            "${right.scope ?: ''}\u0000${right.id ?: ''}\u0000${right.kind ?: ''}\u0000${right.path ?: ''}"
    }
    def cohortOutputs = new TreeMap(recordList
        .findAll { record -> record.scope == 'cohort' }
        .groupBy { record -> record.kind }
        .collectEntries { kind, outputRecords ->
            [kind, outputRecords.collect { record -> record.path }.findAll { path -> path }.unique().sort()]
        })
    def output = [
        schema_version: '2.0',
        pipeline: 'imprint',
        pipeline_version: manifest_info.pipeline_version,
        cohort: manifest_info.cohort,
        mode: manifest_info.mode,
        library_mode: manifest_info.library_mode,
        dupcaller_version: manifest_info.dupcaller_version,
        barcode_chemistry: manifest_info.barcode_chemistry,
        samples: groupManifestRecords(recordList, 'sample'),
        pairs: groupManifestRecords(recordList, 'pair'),
        donors: groupManifestRecords(recordList, 'donor'),
        cohort_outputs: cohortOutputs
    ]
    def json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(jsonSafe(output)))
    task.workDir.resolve('output_manifest.json').text = json + '\n'
}

process EMIT_PROVENANCE {
    label 'process_single'
    tag 'cohort'
    publishDir "${params.outdir}/cohort", mode: 'copy'
    executor 'local'

    input:
    val run_info

    output:
    path 'run_params.json', emit: provenance

    exec:
    def output = [schema_version: '2.0', pipeline: 'imprint'] + run_info
    def json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(jsonSafe(output)))
    task.workDir.resolve('run_params.json').text = json + '\n'
}
