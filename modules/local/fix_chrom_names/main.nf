process FIX_CHROM_NAMES_EXTRACT_PASS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::htslib=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data' :
        'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.chrFixed.vcf.gz"),                                                                                                 emit: vcf
    tuple val(meta), path("*.chrFixed.vcf.gz.tbi"),                                                                                             emit: tbi
    tuple val("${task.process}"), val('perl'), eval("perl --version | grep 'version' | sed 's/.*(v//g' | sed 's/).*//g'"),  topic: versions,    emit: versions_perl
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"),                  topic: versions,    emit: versions_tabix
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"),                                topic: versions,    emit: versions_bgzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fix_chrom_names_extract_pass.sh ${vcf} ${prefix}.chrFixed.vcf

    bgzip -c ${prefix}.chrFixed.vcf > ${prefix}.chrFixed.vcf.gz
    tabix -p vcf ${prefix}.chrFixed.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.chrFixed.vcf.gz 
    touch ${prefix}.chrFixed.vcf.gz.tbi
    """

}