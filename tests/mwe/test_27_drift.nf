// MWE for issue #27 — confirm PHARMCAT_VCF_PROCESSING joins channels by
// `meta.id`, even when the two channels carry meta dicts that share `id`
// but differ in other fields.
//
// At runtime the VCF channel (post PHARMCAT_VCFPREPROCESSOR) and the BED
// channel (from QC_BAM) come from different upstream subworkflows that
// mutate meta differently. A `.join` on the *full* meta map therefore
// fails to match; with `.failOnMismatch:true` it errors, with `false` it
// silently drops or mispairs samples.
//
// This MWE forces the condition deterministically:
//   - VCF channel meta carries an extra `data_type` field
//   - BED channel meta has only `id`
//   - BED emission order is REVERSED — would fall back to positional
//     pairing under the pre-fix `.map { meta, bed -> [bed] }` pattern
//
// The runner script (tests/mwe/run_demo.sh) asserts:
//   POST-FIX: sampleA's BCFTOOLS_VIEW gets A_target.pass.bed (and B gets B's)
//   PRE-FIX:  silent mispair — sampleA gets B_target.pass.bed, vice versa.

include { PHARMCAT_VCF_PROCESSING } from '../../subworkflows/local/pharmcat_vcf_processing'

workflow {
    ch_vcf = Channel.of(
        [[id:'sampleA', data_type:'fastq_gz'],
         file("${projectDir}/inputs/sampleA.vcf.gz"),
         file("${projectDir}/inputs/sampleA.vcf.gz.tbi")],
        [[id:'sampleB', data_type:'fastq_gz'],
         file("${projectDir}/inputs/sampleB.vcf.gz"),
         file("${projectDir}/inputs/sampleB.vcf.gz.tbi")]
    )

    // BED channel: meta has only `id` (no data_type — drift!), and emission
    // order is B then A to expose positional pairing of the pre-fix code.
    ch_bed = Channel.of(
        [[id:'sampleB'], file("${projectDir}/inputs/B_target.pass.bed")],
        [[id:'sampleA'], file("${projectDir}/inputs/A_target.pass.bed")]
    )

    ch_fasta = Channel.value([[id:'ref'], file("${projectDir}/inputs/reference.fna.bgz")])
    ch_idx   = Channel.value([[id:'ref'], file("${projectDir}/inputs/reference.fna.bgz.fai")])
    ch_pos   = Channel.value([
        [id:'pos'],
        file("${projectDir}/inputs/positions.vcf.gz"),
        file("${projectDir}/inputs/positions.vcf.gz.tbi")
    ])
    ch_uni   = Channel.value([
        [id:'uni'],
        file("${projectDir}/inputs/uni.vcf.gz"),
        file("${projectDir}/inputs/uni.vcf.gz.tbi")
    ])

    PHARMCAT_VCF_PROCESSING(ch_vcf, ch_fasta, ch_idx, ch_pos, ch_uni, ch_bed)
}
