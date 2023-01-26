process odgiLayout {
  publishDir "${params.outdir}/odgi_layout", mode: "${params.publish_dir_mode}"
  input:
  path(graph)

  output:
  tuple path(graph), path("${graph}.lay")

  """
  odgi layout \
    -i $graph \
    -o ${graph}.lay \
    -t ${task.cpus} -P
  """
}

process odgiDraw {
  publishDir "${params.outdir}/odgi_draw", mode: "${params.publish_dir_mode}"

  input:
  tuple path(graph), path(layoutGraph)

  output:
  path("${graph}.draw_multiqc.png")

  """
  odgi draw \
    -i $graph \
    -c $layoutGraph \
    -p ${graph}.draw_multiqc.png \
    -C \
    -w 20 \
    -H 1000 -t ${task.cpus}
  odgi draw \
    -i $graph \
    -c $layoutGraph \
    -p ${graph}.draw.png \
    -H 100 -t ${task.cpus}
  """
}

workflow ODGI_2D {
  take:
  graph
  
  main:
  
  odgiLayout(graph)
  odgiDraw(odgiLayout.out)

  emit:
  draw = odgiDraw.out
}
