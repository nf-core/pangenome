/*
** PLEASE READ THIS BEFORE COMPLAINING ABOUT THIS NOT PERFECTLY NF-CORE COMPLIANT MODULE!
**
*** Why did you not use the official VG_DECONSTRUCT module https://github.com/nf-core/modules/tree/master/modules/nf-core/vg/deconstruct?
*** - This official module uses vg 1.43.0, however, due to a bug https://github.com/vgteam/vg/issues/3807 in vg, I would require vg 1.40.0.
*** - Since there is no older module version available on nf-core/modules, I also can't go back to an older commit with the version I need. So here we are.
**
*** Why did you chose this custom docker container and not a container from Bioconda?
*** - The current pggb Bioconda container does not include vcfbub and vcfwave, because both these tools are currently not available on Bioconda.
**
*** Why did you not split up all the commands [vg deconstruct, bcftools stats, bcftools annotate, bgzip, vcfbub, vcfwave] into several modules?
*** - Some of the tools are not available on Bioconda, so it's hard to generate proper modules.
*** - The whole VG_DECONSTRUCT workflow is quite complex and hard to maintain. I basically copied it over from https://github.com/pangenome/pggb/blob/master/pggb#L588-L618.
*** - The original code is still work in progress, so it could change quickly and drastically. Which would mean I have to add or delete or modify at least 7 modules just to get this to run
*** - It is much easier to maintain like this.
**
** I am aware that this module is currently not fitting nf-core standards, but I think I gave good reasons why. Once the code in pggb becomes stable and all tools
** available on Bioconda the goal is to enhance this module so that it actually fits the standards. Thanks!
**
*/
process VG_DECONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

/*  ATTENTION: CURRENTLY VCFBUB AND VCFWAVE ARE NOT AVAILABLE ON BIOCONDA SO I HAVE TO USE THIS IMAGE.

    conda "bioconda::pggb=0.5.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pggb:0.5.3--hdfd78af_2':
        'quay.io/biocontainers/pggb:0.5.3--hdfd78af_2' }"
*/
    container "ghcr.io/pangenome/pggb:20240313103308d2dc38"

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
    if [[ "$vcf_spec" == *":"* ]]; then
        pop_length=\$(echo "$vcf_spec" | cut -f 2 -d:)
    else
        pop_length=""
    fi

    if [[ -z \$pop_length ]]; then
        pop_length=0
    fi

    vcf="${graph}".\$(echo \$ref | tr '/|' '_').vcf
    vg deconstruct -P \$ref -H "#" -e -a -t "${task.cpus}" "${graph}" > \$vcf
    bcftools stats \$vcf > \$vcf.stats

    if [[ \$pop_length -gt 0 ]]; then
        vcf_decomposed=${graph}.final.\$(echo \$ref | tr '/|' '_').decomposed.vcf
        vcf_decomposed_tmp=\$vcf_decomposed.tmp.vcf
        bgzip -c -@ ${task.cpus} \$vcf > \$vcf.gz
        vcfbub -l 0 -a \$pop_length --input \$vcf.gz | vcfwave -I 1000 -t ${task.cpus} > \$vcf_decomposed_tmp
        #TODO: to remove when vcfwave will be bug-free
        # The TYPE info sometimes is wrong/missing
        # There are variants without the ALT allele
        bcftools sort \$vcf_decomposed_tmp | bcftools annotate -x INFO/TYPE | awk '\$5 != "."' > \$vcf_decomposed
        rm \$vcf_decomposed_tmp \$vcf.gz
        bcftools stats \$vcf_decomposed > \$vcf_decomposed.stats
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pggb: \$(pggb --version 2>&1 | grep -o 'pggb .*' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
