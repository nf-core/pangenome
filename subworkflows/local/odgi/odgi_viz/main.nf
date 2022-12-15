process odgiViz {
  publishDir "${params.outdir}/odgi_viz", mode: "${params.publish_dir_mode}"

  input:
    path(graph)

  output:
    path("${graph}.viz*.png")

  script:
    """
    odgi viz -i $graph -o ${graph}.viz_multiqc.png -x 1500 -y 500 -a 10 -I ${params.smoothxg_consensus_prefix}
    odgi viz -i $graph -o ${graph}.viz_pos_multiqc.png -x 1500 -y 500 -a 10 -I ${params.smoothxg_consensus_prefix} -u -d
    odgi viz -i $graph -o ${graph}.viz_depth_multiqc.png -x 1500 -y 500 -a 10 -I ${params.smoothxg_consensus_prefix} -m
    odgi viz -i $graph -o ${graph}.viz_inv_multiqc.png -x 1500 -y 500 -a 10 -I ${params.smoothxg_consensus_prefix} -z
    odgi viz -i $graph -o ${graph}.viz_O_multiqc.png -x 1500 -y 500 -a 10 -I ${params.smoothxg_consensus_prefix} -O
    """
}

workflow ODGI_VIZ {
  take:
  graph
  
  main:
  
  odgiViz(graph)

  emit:
  viz = odgiViz.out
}
