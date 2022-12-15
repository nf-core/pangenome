process odgiStats {
  publishDir "${params.outdir}/odgi_stats", mode: "${params.publish_dir_mode}"

  input:
    path(graph)

  output:
    path("${graph}.stats.yaml")

  """
  odgi stats -i "${graph}" -m --threads ${task.cpus} > "${graph}.stats.yaml" 2>&1
  """
}

workflow ODGI_STATS {
  take:
  graph
  
  main:
  
  odgiStats(graph)

  emit:
  stats = odgiStats.out
}
