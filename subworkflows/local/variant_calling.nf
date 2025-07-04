include { SENTIEON_HAPLOTYPING                      } from '../../modules/local/sentieon/main'
include { GATK_HAPLOTYPING                          } from '../../modules/local/gatk/main'
include { FREEBAYES                                 } from '../../modules/local/freebayes/main'
include { BCFTOOLS_CALL                             } from '../../modules/local/bcftools/main'
include { VT_DECOMPOSE_NORMALIZE as NORM_SENTIEON   } from '../../modules/local/vt/main'
include { VT_DECOMPOSE_NORMALIZE as NORM_GATK       } from '../../modules/local/vt/main'
include { VT_DECOMPOSE_NORMALIZE as NORM_FREEBAYES  } from '../../modules/local/vt/main'
include { VT_DECOMPOSE_NORMALIZE as NORM_BCFTOOLS   } from '../../modules/local/vt/main'
include { AGGREGATE_VCFS                            } from '../../modules/local/aggregate_vcfs/main'



workflow VARIANT_CALLING {

    take:
        bam_ch

    main:
        ch_versions = Channel.empty()

        // Sentieon haplotyping
        SENTIEON_HAPLOTYPING ( bam_ch )
        ch_versions = ch_versions.mix(SENTIEON_HAPLOTYPING.out.versions)

        NORM_SENTIEON ( SENTIEON_HAPLOTYPING.out.haplotypes, "sentieon" ).set { ch_sentieon_vcf }
        ch_versions = ch_versions.mix(ch_sentieon_vcf.versions)
        
        // GATK haplotyping
        GATK_HAPLOTYPING ( bam_ch )
        ch_versions = ch_versions.mix(GATK_HAPLOTYPING.out.versions)

        NORM_GATK ( GATK_HAPLOTYPING.out.haplotypes, "gatk" ).set { ch_gatk_vcf }
        ch_versions = ch_versions.mix(ch_gatk_vcf.versions)

        // Freebayes haplotyping
        FREEBAYES ( bam_ch )
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)

        NORM_FREEBAYES ( FREEBAYES.out.freebayes_vcf, "freebayes" ).set { ch_freebayes_vcf }
        ch_versions = ch_versions.mix(ch_freebayes_vcf.versions)

        // samtools mpileup
        BCFTOOLS_CALL ( bam_ch )
        ch_versions = ch_versions.mix(BCFTOOLS_CALL.out.versions)

        NORM_BCFTOOLS ( BCFTOOLS_CALL.out.bcftools_vcf, "bcftools" ).set { ch_bcftools_vcf }
        ch_versions = ch_versions.mix(ch_bcftools_vcf.versions)


        // Aggregate all callers to one VCF
        ch_vcf = Channel.empty().mix(
                ch_sentieon_vcf.decomposed_normalized_vcfs, 
                ch_gatk_vcf.decomposed_normalized_vcfs, 
                ch_freebayes_vcf.decomposed_normalized_vcfs, 
                ch_bcftools_vcf.decomposed_normalized_vcfs
            ).groupTuple(by: [0,1])

        AGGREGATE_VCFS ( ch_vcf )
        // ch_versions = ch_versions.mix(AGGREGATE_VCFS.out.versions)



    emit:
        sentieon_vcf        =   ch_sentieon_vcf.decomposed_normalized_vcfs
        gatk_vcf            =   ch_gatk_vcf.decomposed_normalized_vcfs
        freebayes_vcf       =   ch_freebayes_vcf.decomposed_normalized_vcfs
        bcftools_vcf        =   ch_bcftools_vcf.decomposed_normalized_vcfs
        aggregate_vcf       =   AGGREGATE_VCFS.out.vcf_agg
        aggregate_vcf_tbi   =   AGGREGATE_VCFS.out.vcf_agg_tbi
        versions            =   ch_versions

}