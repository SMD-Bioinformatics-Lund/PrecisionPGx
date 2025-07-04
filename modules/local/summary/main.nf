process SAMPLE_TARGET_LIST {
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(detected_variants)

    output:
        tuple val(group), val(meta), file("*.pgx_target_interval.list"),    emit: pgx_target_interval_list
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        reform_genomic_region.py \
            --output_file=${prefix}.pgx_target_interval.list \
            --detected_variants=$detected_variants \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pgx_target_interval.list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process DEPTH_OF_TARGETS {
    // Get read depth of variant locations at wildtrype-called positions
    label 'process_low'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai)

    output:
        tuple val(group), val(meta), file("*.pgx_depth_at_missing.gdf"),    emit: pgx_depth_at_missing
        path "versions.yml",                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}"
        def avail_mem   = 6144
        if (!task.memory) {
            log.info '[GATK DepthOfCoverage] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        java -Xmx${avail_mem}M -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage $args \
        -I $bam \
        -o ${prefix}.pgx_depth_at_missing.gdf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(java -jar /usr/GenomeAnalysisTK.jar --version)
        END_VERSIONS
        """

    stub:
        def prefix      = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pgx_depth_at_missing.gdf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(java -jar /usr/GenomeAnalysisTK.jar --version)
        END_VERSIONS
        """
}

process GET_PADDED_BAITS {
    label 'process_single'
    label 'stage'
    tag "Get_padded_baits_list"

    output:
        path "*padded_bait_interval.list",  emit: padded_baits_list
        path "versions.yml",                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def prefix      = task.ext.prefix ?: ''
        """
        reform_genomic_region.py \
            --output_file=${prefix}padded_bait_interval.list \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix      = task.ext.prefix ?: ''
        """
        touch ${prefix}padded_bait_interval.list

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process DEPTH_OF_BAITS {
    // Get read depth of baits
    label 'process_low'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(bam), file(bai) 

    output:
        tuple val(group), val(meta), file("*.pgx.gdf"), emit: padded_baits_list
        path "versions.yml",                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args        = task.ext.args   ?: ''
        def prefix      = task.ext.prefix ?: "${meta.group}"
        def avail_mem   = 6144
        if (!task.memory) {
            log.info '[GATK DepthOfCoverage] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        # NOTE: does not work with openjdk-11, openjdk-8 works
        java -Xmx${avail_mem}M -jar /usr/GenomeAnalysisTK.jar -T DepthOfCoverage $args \
        -I $bam \
        -o ${prefix}.pgx.gdf 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(java -jar /usr/GenomeAnalysisTK.jar --version)
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pgx.gdf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(java -jar /usr/GenomeAnalysisTK.jar --version)
        END_VERSIONS
        """
}

process PADDED_BED_INTERVALS {
    label 'process_single'
    label 'stage'
    tag "pgx_bed_padded_intervals"

    output:
        path '*padded_bait_interval.bed',   emit: padded_bed_intervals
        path "versions.yml",                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: ''
        """
        reform_genomic_region.py \
            --output_file=${prefix}padded_bait_interval.bed \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: ''
        """
        touch ${prefix}padded_bait_interval.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}

process APPEND_ID_TO_GDF {
    //  Add variant id to appropriate location in gdf
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(gdf)

    output:
        tuple val(group), val(meta), file("*.pgx_depth_at_missing_annotated.gdf"),  emit: depth_at_missing_annotate_gdf
        path "versions.yml",                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        append_rsid_to_gdf.py \
            --input_gdf=$gdf \
            --output_file=${prefix}.pgx_depth_at_missing_annotated.gdf \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pgx_depth_at_missing_annotated.gdf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed -e 's/Python //g')
        END_VERSIONS
        """
}


process EXTRACT_CYP2D6_HAPLOTYPES {
    label 'process_low'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(match_json)

    output:
        tuple val(group), val(meta), file("*.cyp2d6.genotypes.txt"),    emit: cy2d6_genotypes

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        grep -wC3 "CYP2D6"  $match_json | tail -n1 | cut -f2 -d ":" | tr -d " " | tr -d "," > ${prefix}.cyp2d6.genotypes.txt
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.cyp2d6.genotypes.txt
        """

}