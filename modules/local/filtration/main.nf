process VARIANT_FILTRATION {
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf), file(tbi)

    output:
        tuple val(group), val(meta), file("*.filtered.haplotypes.vcf.gz"), file("*.filtered.haplotypes.vcf.gz.tbi"),    emit: haplotypes_filtered
        path "versions.yml",                                                                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        gunzip -c $vcf > ${prefix}.haplotypes.vcf
        variant_filtration.py \
            --input_vcf=${prefix}.haplotypes.vcf \
            $args \
            --output_file=${prefix}.filtered.haplotypes.vcf
        
        bgzip -c ${prefix}.filtered.haplotypes.vcf > ${prefix}.filtered.haplotypes.vcf.gz
        tabix -p vcf ${prefix}.filtered.haplotypes.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.filtered.haplotypes.vcf.gz ${prefix}.filtered.haplotypes.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
            bgzip: \$(bgzip --v | grep 'bgzip' | sed 's/.* //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """

}