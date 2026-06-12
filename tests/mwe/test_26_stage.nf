// MWE for issue #26 — confirm `stageInMode = 'copy'` is in effect for
// PHARMCAT_VCFPREPROCESSOR, so htslib's touch/rewrite of
// `reference.fna.bgz.fai` cannot reach the shared source file via a symlink.
//
// Invokes the subworkflow (not the module directly) so the project's
// `withName: '.*PHARMCAT_VCF_PROCESSING:PHARMCAT_VCFPREPROCESSOR'`
// directive — which carries the fix — actually applies.
//
// The runner script (tests/mwe/run_demo.sh) asserts:
//   POST-FIX: PHARMCAT_VCFPREPROCESSOR workdir contains 0 input symlinks,
//             every staged file is a real per-task copy.
//   PRE-FIX:  every staged input is a symlink back to projectDir/inputs/,
//             so any htslib write through `reference.fna.bgz.fai` modifies
//             the shared source — the race condition the issue describes.

include { PHARMCAT_VCF_PROCESSING } from '../../subworkflows/local/pharmcat_vcf_processing'

workflow {
    ch_vcf = Channel.of([
        [id:'sampleA'],
        file("${projectDir}/inputs/sampleA.vcf.gz"),
        file("${projectDir}/inputs/sampleA.vcf.gz.tbi")
    ])
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
    ch_bed   = Channel.of([[id:'sampleA'], file("${projectDir}/inputs/A_target.pass.bed")])

    PHARMCAT_VCF_PROCESSING(ch_vcf, ch_fasta, ch_idx, ch_pos, ch_uni, ch_bed)
}
