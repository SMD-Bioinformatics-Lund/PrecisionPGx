process VT_DECOMPOSE_NORMALIZE {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), file(vcf), file(tbi)
        val(vc)

    output:
        tuple val(group), val(meta), val(vc), file("*_${vc}.vcf.gz"), file("*_${vc}.vcf.gz.tbi"),   emit: decomposed_normalized_vcfs
        path "versions.yml",                                                                        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${group}"
        """
        vt decompose ${vcf} -o ${prefix}.${vc}.decomposed.vcf.gz
        vt normalize  ${prefix}.${vc}.decomposed.vcf.gz $args | vt uniq - -o ${prefix}_${vc}.vcf.gz
        tabix -p vcf ${prefix}_${vc}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vt-decompose: \$(echo \$(vt decompose 2>&1) | sed 's/.*decompose v//; s/ .*//')
            vt-normalize: \$(echo \$(vt normalize 2>&1) | sed 's/.*normalize v//; s/ .*//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        touch ${prefix}_${vc}.vcf.gz  ${prefix}_${vc}.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            vt-decompose: \$(echo \$(vt decompose 2>&1) | sed 's/.*decompose v//; s/ .*//')
            vt-normalize: \$(echo \$(vt normalize 2>&1) | sed 's/.*normalize v//; s/ .*//')
        END_VERSIONS
        """
}