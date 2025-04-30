process SPLIT_APPROX_MAPPINGS_IN_CHUNKS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::wfmash=0.12.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wfmash:0.12.6--h11f254b_1':
        'biocontainers/wfmash:0.12.6--h11f254b_1' }"

    input:
    tuple val(meta), path(paf)

    output:
    tuple val(meta), path("*.paf"), emit: pafs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    split_approx_mappings_in_chunks.py $paf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pggb: \$(pggb --version 2>&1 | grep -o 'pggb .*' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
