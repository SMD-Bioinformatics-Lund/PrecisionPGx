process GATK_HAPLOTYPING {
    label 'process_medium_cpus'
    label 'process_medium_memory'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.haplotypes.vcf.gz"), file("*.haplotypes.vcf.gz.tbi"),  emit: haplotypes
        path "versions.yml",                                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}.GATK"
        def avail_mem   = 12288
        if (!task.memory) {
            log.info '[GATK CollectAllelicCounts] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        gatk --java-options "-Xmx${avail_mem}M" HaplotypeCaller $args -I $bam -O ${prefix}.haplotypes.vcf
        bgzip -c ${prefix}.haplotypes.vcf > ${prefix}.haplotypes.vcf.gz
        tabix ${prefix}.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix      = task.ext.prefix ?: "${meta.group}.GATK"
        """
        touch ${prefix}.haplotypes.vcf.gz ${prefix}.haplotypes.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
}
