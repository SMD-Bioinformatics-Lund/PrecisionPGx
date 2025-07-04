process GET_CLIINICAL_GUIDELINES {
    // Given detected variants, get possible Haplotype combinations
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(detected_variants)

    output:
        tuple val(group), val(meta), file("*.possible_diplotypes.tsv"), emit: possible_diplotypes
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        get_possible_diplotypes.py \
            --variant_csv $detected_variants \
            --output ${prefix}.possible_diplotypes.tsv  \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.possible_diplotypes.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process GET_INTERACTION_GUIDELINES {
    // Given Haplotype Combinations, get possible interactions betweens these
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(possible_diplotypes)

    output:
        tuple val(group), val(meta), file("*.possible_interactions.tsv"),   emit: possible_interactions
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        get_interaction_guidelines.py \
            --diploids $possible_diplotypes \
            --output ${prefix}.possible_interactions.tsv \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.possible_interactions.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}