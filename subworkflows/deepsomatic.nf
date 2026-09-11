include { CALL } from '../modules/deepsomatic'

workflow DEEPSOMATIC {
    take:
        ch_paired_bams // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

    main:
        CALL(ch_paired_bams)

    emit:
        CALL.out.vcf
}
