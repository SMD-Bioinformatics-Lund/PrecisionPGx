process DETECTED_VARIANTS {
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf) 
        
    output:
        tuple val(group), val(meta), file("*.detected_variants.tsv"),   emit: detected_tsv
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        get_target_variants.py \
            --vcf $vcf \
            --output ${prefix}.detected_variants.tsv \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.detected_variants.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}