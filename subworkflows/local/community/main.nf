def wfmash_prefix = "wfmash.map.community"

include { WFMASH_MAP } from '../wfmash_map/main'

workflow COMMUNITY {
  take:
  ch_fasta
  fai_path
  gzi_path

  main:

  WFMASH_MAP(ch_fasta, fai_path, gzi_path, wfmash_prefix)

  emit:
  ch_fasta

}