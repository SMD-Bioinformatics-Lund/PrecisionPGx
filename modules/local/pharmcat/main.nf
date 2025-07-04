process PHARMCAT_PREPROCESSING {
    label 'process_single'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf), file(tbi)

    output:
        tuple val(group), val(meta), file("*.pharmcat.preprocessed.vcf"),   emit: pharmcat_preprocessed_vcf
        // tuple val(group), val(meta), file("*.pharmcat.pgx_regions.vcf.bgz"), file("*.pharmcat.pgx_regions.vcf.bgz.csi"),    emit: pharmcat_pgx_regions
        // tuple val(group), val(meta), file("*.pharmcat.normalized.vcf.bgz"), file("*.pharmcat.normalized.vcf.bgz.csi"),      emit: pharmcat_normalized
        path "versions.yml",                                                                                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        pharmcat_vcf_preprocessor.py -vcf $vcf --base-filename ${prefix}.pharmcat $args
        gunzip -c ${prefix}.pharmcat.preprocessed.vcf.bgz > ${prefix}.pharmcat.preprocessed.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pharmcat_vcf_preprocessor: \$(pharmcat_vcf_preprocessor.py -V 2>&1 | sed -e 's/PharmCAT VCF Preprocessor //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pharmcat.preprocessed.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pharmcat_vcf_preprocessor: \$(pharmcat_vcf_preprocessor.py -V 2>&1 | sed -e 's/PharmCAT VCF Preprocessor //g')
        END_VERSIONS
        """

}


process PHARMCAT_RUN {
    label 'process_very_few_cpus'
    label 'process_medium_memory'
    label 'stage'
    tag "$meta.group"

    input:
        tuple val(group), val(meta), file(vcf)

    output:
        tuple val(group), val(meta), file("*.phenotype.json"),  emit: pharmcat_pheno_json
        tuple val(group), val(meta), file("*.match.json"),      emit: pharmcat_match_json
        tuple val(group), val(meta), file("*.report.html"),     emit: pharmcat_report,              optional: true
        tuple val(group), val(meta), file("*.report.json"),     emit: pharmcat_report_json,         optional: true
        tuple val(group), val(meta), file("*.tsv"),             emit: pharmcat_calls_tsv,           optional: true  
        path "versions.yml",                                    emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args    = task.ext.args   ?: ''
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        pharmcat -vcf  $vcf --base-filename ${prefix}.pharmcat $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            PharmCAT: \$(pharmcat --version 2>&1 | sed -e 's/PharmCAT //g')
        END_VERSIONS
        """

    stub:
        def prefix  = task.ext.prefix ?: "${meta.group}"
        """
        touch ${prefix}.pharmcat.phenotype.json
        touch ${prefix}.match.json
        touch ${prefix}.report.html
        touch ${prefix}.report.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            PharmCAT: \$(pharmcat --version 2>&1 | sed -e 's/PharmCAT //g')
        END_VERSIONS
        """

}


// process PHARMCAT_VCF_BED {
//     label 'process_medium'
//     label 'stage'
//     tag "pharmcat_vcf_bed"

//     output:
//         path "*.bed",               emit: pharmcat_positions_bed
//         path "versions.yml",        emit: versions

//     when:
//         task.ext.when == null || task.ext.when

//     script:
//         def args    = task.ext.args   ?: ''
//         def prefix  = task.ext.prefix ?: "pharmcat_vcf_positions"
//         """
//         create_pharmcat_bed.py --output ${prefix}.bed $args

//         cat <<-END_VERSIONS > versions.yml
//         "${task.process}":
//             python: \$(python --version | sed 's/.*\s//g')
//         END_VERSIONS
//         """

//     stub:
//         def prefix  = task.ext.prefix ?: "pharmcat_vcf_positions"
//         """
//         create_pharmcat_bed.py --output ${prefix}.bed $args

//         cat <<-END_VERSIONS > versions.yml
//         "${task.process}":
//             python: \$(python --version | sed 's/.*\s//g')
//         END_VERSIONS
//         """

// }