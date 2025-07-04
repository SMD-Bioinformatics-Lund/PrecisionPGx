include { VCFTOOLS_FILTER                      } from '../../modules/local/vcftools/main'



workflow VARIANT_FILTRATION {

    take:
        vcf

    main:
        ch_versions = Channel.empty()

        // Sentieon haplotyping
        VCFTOOLS_FILTER ( vcf )
        ch_versions = ch_versions.mix(VCFTOOLS_FILTER.out.versions)


    emit:
        filtered_vcf    =   VCFTOOLS_FILTER.out.filtered_vcf
        versions        =   ch_versions

}