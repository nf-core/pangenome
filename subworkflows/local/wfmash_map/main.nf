// We can't change global parameters inside this scope, so we build the ones we need locally
def wfmash_merge_cmd = params.wfmash_merge_segments ? "-M" : ""
def wfmash_exclude_cmd = params.wfmash_exclude_delim ? "-Y ${params.wfmash_exclude_delim}" : "-X"
def wfmash_split_cmd = params.wfmash_no_splits ? "-N" : ""
def wfmash_block_length = params.wfmash_segment_length*5
def wfmash_block_length_cmd = "-l ${wfmash_block_length}"
def wfmash_mash_kmer_cmd = "-k ${params.wfmash_mash_kmer}"
def wfmash_kmer_thres_cmd = "-H ${params.wfmash_mash_kmer_thres}"
def wfmash_n_mappings_minus_1 = params.n_haplotypes - 1
def wfmash_sparse_map_cmd = ""
if (params.wfmash_sparse_map == "auto") {
  n = n_haps
  x = Math.log(n)/n * 10
  wfmash_sparse_map_frac = 1
  if (x >= 1) {
    wfmash_sparse_map_frac = x
  }
  wfmash_sparse_map_cmd = "-x${wfmash_sparse_map_frac}"
} else {
  if (params.wfmash_sparse_map != null) {
    wfmash_sparse_map_cmd = "-x${params.wfmash_sparse_map}"
  }
}
def wfmash_temp_dir = params.wfmash_temp_dir ? "-B${params.wfmash_temp_dir}" : ""

def make_file_prefix = { f -> """\
${f.getName()}\
""" }

process wfmashMap {
  publishDir "${params.outdir}/wfmash_map", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta)
    path(fai)
    path(gzi)
    val(wfmash_prefix)

  output:
    tuple val(f), path("${f}.${wfmash_prefix}.map.paf")

  """
  wfmash ${wfmash_exclude_cmd} \
     -s ${params.wfmash_segment_length} \
     ${wfmash_block_length_cmd} \
     ${wfmash_merge_cmd} \
     ${wfmash_split_cmd} \
     ${wfmash_mash_kmer_cmd} \
     ${wfmash_kmer_thres_cmd} \
     ${wfmash_sparse_map_cmd} \
     -p ${params.wfmash_map_pct_id} \
     -n ${wfmash_n_mappings_minus_1} \
     ${wfmash_temp_dir} \
     -t ${task.cpus} \
     -m \
     $fasta $fasta \
     >${f}.${wfmash_prefix}.map.paf
  """  
}

workflow WFMASH_MAP {
  take:
  ch_fasta
  fai_path
  gzi_path
  prefix

  main:

  fasta = ch_fasta.map { f -> tuple(make_file_prefix(f), f) }
  ch_wfmash_map = wfmashMap(fasta, fai_path, gzi_path, prefix)

  emit:
  ch_wfmash_map
}