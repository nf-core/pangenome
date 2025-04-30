process NET2COMMUNITIES {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pggb=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pggb:0.5.4--hdfd78af_0':
        'biocontainers/pggb:0.5.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(txts)

    output:
    path("*.community.*.txt"), emit: communities
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    net2communities.py \
    -e ${prefix}.paf.edges.list.txt \
    -w ${prefix}.paf.edges.weights.txt \
    -n ${prefix}.paf.vertices.id2name.txt \
    --accurate-detection \
    --output-prefix ${prefix} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pggb: \$(pggb --version 2>&1 | grep -o 'pggb .*' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
