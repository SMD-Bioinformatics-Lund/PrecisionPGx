process ONTARGET_BAM {
    label 'process_low'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group),  val(meta), file("*.dedup.ontarget.pgx.bam"), file("*.dedup.ontarget.pgx.bam.bai"),   emit: bam_ontarget
        path "versions.yml",                                                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        samtools view $args -b $bam > ${prefix}.dedup.ontarget.pgx.bam
        samtools index ${prefix}.dedup.ontarget.pgx.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | grep 'samtools' | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.dedup.ontarget.pgx.bam ${prefix}.dedup.ontarget.pgx.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | grep 'samtools' | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process ONTARGET_VCF {
    label 'process_low'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf), file(tbi) 

    output:
        tuple val(group), val(meta), file("*.ontarget.filtered.haplotypes.vcf.gz"), file("*.ontarget.filtered.haplotypes.vcf.gz.tbi"),  emit: vcf_ontarget
        path "versions.yml",                                                                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        bcftools view $args -o ${prefix}.ontarget.filtered.haplotypes.vcf $vcf
        bgzip -c ${prefix}.ontarget.filtered.haplotypes.vcf > ${prefix}.ontarget.filtered.haplotypes.vcf.gz
        tabix -p vcf ${prefix}.ontarget.filtered.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | grep 'bcftools' 2>&1 | sed 's/^.*bcftools //')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.ontarget.filtered.haplotypes.vcf.gz
        touch ${prefix}.ontarget.filtered.haplotypes.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version | grep 'bcftools' 2>&1 | sed 's/^.*bcftools //')
        END_VERSIONS
        """
}


