ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)

def make_file_prefix = { f -> """\
${f.getName()}\
""" }

process odgiSqueeze {
  publishDir "${params.outdir}/odgi_squeeze", mode: "${params.publish_dir_mode}"

  input:
    path(graphs)
    tuple val(f), path(fasta)

  output:
    path("${f}.squeezed.og"), emit: squeezed

  """
  ls *.og > files
  odgi squeeze -f files --optimize -o "${f}".squeezed.og --threads ${task.cpus}
  """
}

process odgi2Gfa {
  publishDir "${params.outdir}/odgi2gfa", mode: "${params.publish_dir_mode}"

  input:
    path(graph)

  output:
    path("*.og.gfa"), emit: gfa

  """
  odgi view -i ${graph} -g > ${graph}.gfa
  """
}

process multiQC {
  publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

  input:
  path odgi_stats
  path odgi_viz
  path odgi_draw
  path(multiqc_config)

  output:
  path "*multiqc_report.html", emit: report
  path "*_data"              , emit: data
  path "*_plots"             , optional:true, emit: plots

  """
  multiqc -s . -c ${multiqc_config}
  """
}

include { ODGI_STATS } from '../odgi/odgi_stats/main'
include { ODGI_VIZ } from '../odgi/odgi_viz/main'
include { ODGI_2D } from '../odgi/odgi2d/main'

workflow SQUEEZE {
  take:
  ch_pggb
  ch_fasta

  main:
  fasta = ch_fasta.map { f -> tuple(make_file_prefix(f), f) }
  ch_squeezed = odgiSqueeze(ch_pggb, fasta)
  // TODO GFA as optional additional output
  if (params.squeeze_gfa) {
    odgi2Gfa(odgiSqueeze.out)
  }
  ODGI_STATS(odgiSqueeze.out)
  ODGI_VIZ(odgiSqueeze.out)
  ODGI_2D(odgiSqueeze.out)
  multiQC(
    ODGI_STATS.out,
    ODGI_VIZ.out,
    ODGI_2D.out,
    ch_multiqc_config
    )

  emit:
  ch_squeezed
}
