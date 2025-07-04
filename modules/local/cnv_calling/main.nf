process CALCULATE_CNV_CYP2D6  {
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(depthfile)

    output:
        tuple val(group), val(meta), file("*.CNV.txt"),                         emit: cypd26_cnv
        path "versions.yml",                                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        postprocess_cnvs.py -s ${prefix} -d ${depthfile} ${args} -o ${prefix}.CNV.txt 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        python --version
        python3 --version
        which python3
        touch ${prefix}.CNV.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

}


process REPORT_CNV_CYP2D6 {
    label 'process_medium'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(cnvfile)

    output:
        tuple val(group), val(meta), file("*report_filtered.txt"), file ("*CYP2D6_vs_CYP2D7_table.txt"), file("*.png"), emit: cnvreport_cypd26

        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """

        generate_CYP2D6_CNV_report_v2.py -s ${prefix} -c ${cnvfile} $args  -o ./${prefix}
        cp  ./${prefix}/${prefix}_CYP2D6_report/* .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

    stub:

        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        python --version
        which python
        touch ${prefix}.report_filtered.txt
        touch ${prefix}.CYP2D6_vs_CYP2D7_table.txt
        touch ${prefix}..png

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1| sed -e 's/Python //g')
        END_VERSIONS
        """

}

// process CALCULATE_CNV_CYP2D6_PON {

// }