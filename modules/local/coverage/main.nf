process COVERAGE_REPORTS {
    tag "${meta.group}"

    input:
        tuple val(group), val(meta), file(depth_file)

    output:
        tuple val(group), val(meta), file("*.coverage_stats.txt"),      emit: cov_stats
        tuple val(group), val(meta), file("*.depth_annotated.csv"),     emit: cov_annotated
        tuple val(group), val(meta), file("*.coverage_report.html"),    emit: cov_html
        tuple val(group), val(meta), file("*.png"),                     emit: cov_plots
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        coverage_reports_and_plots.py --sample_id ${meta.id} --depth_file ${depth_file} --out_prefix ${prefix} $args 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/.*\s//g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}"
        """
        which python
        python3 --version
        python --version
        touch ${prefix}.coverage_stats.txt
        touch ${prefix}.depth_annotated.csv
        touch ${prefix}.coverage_report.html
        touch ${prefix}.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/.*\s//g')
        END_VERSIONS
        """
}

process PROBLAMATIC_REGIONS_COVERAGE_REPORTS {
    tag "${meta.group}"

    input:
        tuple val(group), val(meta), file(depth_file)

    output:
        tuple val(group), val(meta), file("*.overlap_report.txt"),      emit: overlap_cov_txt
        tuple val(group), val(meta), file("*.overlap_report.html"),     emit: overlap_cov_html
        tuple val(group), val(meta), file("*.png"),                     emit: overlap_cov_plots
        path "versions.yml",                                            emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        problamatic_coverage_calculation.py --sample_name ${meta.id} --depths_file ${depth_file} $args 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/.*\s//g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.overlap_report.txt
        touch ${prefix}.overlap_report.html
        touch ${prefix}.png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/.*\s//g')
        END_VERSIONS
        """
}


process UNCOVERED_PHARMCAT_POSTIONS {
    tag "${meta.group}"

    input:
        tuple val(group), val(meta), file(depth_file)

    output:
        tuple val(group), val(meta), file("*.pharmcat.pos.low.coverage.txt"),   emit: lowcovered_pc_positions
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        pharmcat_coverage_calc.py --sample ${meta.id} --depths_file ${depth_file} $args 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/.*\s//g')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pharmcat.pos.low.coverage.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/.*\s//g')
        END_VERSIONS
        """
}
