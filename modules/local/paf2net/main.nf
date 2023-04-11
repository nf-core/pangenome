process PAF2NET {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pggb=0.5.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pggb:0.5.3--hdfd78af_2':
        'quay.io/biocontainers/pggb:0.5.3--hdfd78af_2' }"

    input:
    tuple val(meta), path(paf)

    output:
    tuple val(meta), path("*.txt"), emit: txts
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    paf2net.py -p $paf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pggb: \$(echo \" 0.5.3\")
    END_VERSIONS
    """
}