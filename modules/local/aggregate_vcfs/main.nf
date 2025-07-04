process AGGREGATE_VCFS {
    label "process_single"
    tag "$group"

    input:
        tuple val(group), val(meta), val(vc), file(vcfs), file(tbis)

    output:
        tuple val(group), val(meta), file("*.agg.vcf"),                                 emit: vcf_agg
        tuple val(group), val(meta), file("*.agg.vcf.gz"), file("*.agg.vcf.gz.tbi"),    emit: vcf_agg_tbi

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${group}"
        """
        vcf_list=\$(echo $vcfs| tr -s " " ",")
        aggregate_vcf.py --vcf \${vcf_list} --sample-order ${meta.id} | vcf-sort -c > ${prefix}.agg.vcf
        bgzip -c ${prefix}.agg.vcf > ${prefix}.agg.vcf.gz
        tabix -p vcf ${prefix}.agg.vcf.gz
        """

    stub:
        def prefix = task.ext.prefix ?: "${group}"
        """
        vcf_list=\$(echo $vcfs| tr -s " " ",")
        echo "vcf_list: \${vcf_list}"
        touch ${group}.agg.vcf
        touch ${group}.agg.vcf.gz
        touch ${group}.agg.vcf.gz.tbi
        """
}
