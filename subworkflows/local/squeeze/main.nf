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
  odgi squeeze -f files --optimize -o "${f}".squeezed.og
  """
}

workflow SQUEEZE {
  take:
  ch_pggb
  ch_fasta

  main:
  fasta = ch_fasta.map { f -> tuple(make_file_prefix(f), f) }
  ch_squeezed = odgiSqueeze(ch_pggb, fasta)
  // TODO GFA as optional additional output
  // TODO stats
  // TODO viz
  // TODO layout
  // TODO MultiQC

  emit:
  ch_squeezed
}
