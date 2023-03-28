process VG_DECONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

/*  ATTENTION: CURRENTLY VCFBUB AND VCFWAVE ARE NOT AVAILABLE ON BIOCONDA SO I HAVE TO USE THIS IMAGE.

    conda "bioconda::pggb=0.5.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pggb:0.5.3--hdfd78af_2':
        'quay.io/biocontainers/pggb:0.5.3--hdfd78af_2' }"
*/
    container "ghcr.io/pangenome/pggb:202303241438332f46d5"

    input:
    tuple val(meta), path(graph), val(vcf_spec)

    output:
    tuple val(meta), path("*.vcf"),                emit: vcf
    path("*.stats"),               optional: true, emit: stats
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ref=\$(echo "$vcf_spec" | cut -f 1 -d:)
    delim=\$(echo "$vcf_spec" | cut -f 2 -d:)
    pop_length=\$(echo "$vcf_spec" | cut -f 3 -d:)

    if [[ -z \$pop_length ]]; then
        pop_length=0
    fi

    vcf="${graph}".\$(echo \$ref | tr '/|' '_').vcf
    vg deconstruct -P \$ref -H \$delim -e -a -t "${task.cpus}" "${graph}" > \$vcf
    bcftools stats \$vcf > \$vcf.stats
    if [[ \$pop_length -gt 0 ]]; then
        vcf_decomposed=${graph}.final.\$(echo \$ref | tr '/|' '_').decomposed.vcf
        vcf_decomposed_tmp=\$vcf_decomposed.tmp.vcf
        bgzip -c -@ ${task.cpus} \$vcf > \$vcf.gz
        vcfbub -l 0 -a \$pop_length --input \$vcf.gz | vcfwave -I 1000 -t ${task.cpus} > \$vcf_decomposed_tmp
        #TODO: to remove when vcfwave will be bug-free
        # The TYPE info sometimes is wrong/missing
        # There are variants without the ALT allele
        bcftools annotate -x INFO/TYPE \$vcf_decomposed_tmp  | awk '\$5 != "."' > \$vcf_decomposed
        rm \$vcf_decomposed_tmp \$vcf.gz
        bcftools stats \$vcf_decomposed > \$vcf_decomposed.stats
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \" 1.40.0\")
    END_VERSIONS
    """
}
