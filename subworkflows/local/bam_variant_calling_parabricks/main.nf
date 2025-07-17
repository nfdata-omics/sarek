//
// DEEPVARIANT germline calling with Parabricks
//

include { PARABRICKS_DEEPVARIANT } from '../../../modules/nf-core/parabricks/deepvariant/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_VARIANT_CALLING_PARABRICKS {
    take:
    cram_intervals // channel: [meta, cram, crai, intervals]
    fasta         // channel: [ref_fasta, fasta]
    dict          // channel: [optional] [ meta, dict ]

    main:
    versions = Channel.empty()

    // Run the Parabricks DeepVariant process
    PARABRICKS_DEEPVARIANT(cram_intervals, fasta)

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_out = PARABRICKS_DEEPVARIANT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more gvcf(s) from the same sample
    gvcf_out = PARABRICKS_DEEPVARIANT.out.gvcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    gvcf_to_merge = gvcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_to_merge = vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    MERGE_DEEPVARIANT_GVCF(gvcf_to_merge, dict)
    MERGE_DEEPVARIANT_VCF(vcf_to_merge, dict)

    // Mix intervals and no_intervals channels together
    gvcf = Channel.empty().mix(MERGE_DEEPVARIANT_GVCF.out.vcf, gvcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.vcf, vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    versions = versions.mix(PARABRICKS_DEEPVARIANT.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_GVCF.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)

    emit:
    vcf
    gvcf

    versions

}