def wfmash_prefix = "wfmash.map.community"
def net2Community_prefix = "${wfmash_prefix}.map.paf"

def make_file_prefix = { f -> """\
${f.getName()}\
""" }

include { WFMASH_MAP } from '../wfmash_map/main'

process paf2Net {
  publishDir "${params.outdir}/paf2net", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(paf)
  output:
    tuple val(f), path("${f}*.txt")
  """
  python3 /paf2net.py -p $paf
  """
}

process net2Communities {
  publishDir "${params.outdir}/net2communities", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(txts)
  output:
    path("${f}.community.*.txt")
  """
  python3 /net2communities.py \
    -e ${f}.${net2Community_prefix}.edges.list.txt \
    -w ${f}.${net2Community_prefix}.edges.weights.txt \
    -n ${f}.${net2Community_prefix}.vertices.id2name.txt \
    --accurate-detection \
    --output-prefix ${f} \
  """
} 

process extractCommunities {
  publishDir "${params.outdir}/extract_communities", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta), path(community)
  output:
    path("${f}*.community.*.fa")
  """
  samtools faidx ${fasta} \$(cat ${community}) > ${community}.fa
  """
}

process bgzip {
  publishDir "${params.outdir}/bgzip", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta), path(extracted_community)
  output:
    path("${f}*.community.*.fa.gz")
  """
  bgzip -@ ${task.cpus} -c ${extracted_community} > ${extracted_community}.gz
  """
}

process samtoolsFaidx {
  publishDir "${params.outdir}/samtools_faidx", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta)

  output:
    path("${f}.fai"), emit: samtools_fai
    path("${f}.gzi"), emit: samtools_gzi

  """
  samtools faidx $fasta
  """
}

workflow COMMUNITY {
  take:
  ch_fasta
  fai_path
  gzi_path

  main:

  if (!fai_path.exists() || !gzi_path.exists()) { // the assumption is that none of the files exist if only one does not exist
    samtoolsFaidx(fasta)
    fai = samtoolsFaidx.out.samtools_fai.collect()
    gzi = samtoolsFaidx.out.samtools_gzi.collect()
  } else {
    fai = channel.fromPath("${params.input}.fai").collect()
    gzi = channel.fromPath("${params.input}.gzi").collect()
  }

  WFMASH_MAP(ch_fasta, fai_path, gzi_path, wfmash_prefix)
  paf2Net(WFMASH_MAP.out)
  net2Communities(paf2Net.out)
  fasta = ch_fasta.map { f -> tuple(make_file_prefix(f), f) }
  extractCommunities(fasta.combine(net2Communities.out.flatten()))
  ch_bgzip_extract_communities = bgzip(fasta.combine(extractCommunities.out.flatten()))

  emit:
  ch_bgzip_extract_communities /// TODO!

}