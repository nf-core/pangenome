ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()

// We can't change global parameters inside this scope, so we build the ones we need locally
def n_haps = 0
if (!params.smoothxg_haplotypes_smooth) {
  n_haps = params.n_haplotypes
}

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

def seqwish_temp_dir = params.seqwish_temp_dir ? "--temp-dir${params.seqwish_temp_dir}" : ""

def smoothxg_block_id_min = params.wfmash_map_pct_id / 100.0
def smoothxg_poa_params_cmd = ""
if (params.smoothxg_poa_params == null) {
  smoothxg_poa_params = "-P 1,19,39,3,81,1"
} else {
  if (params.smoothxg_poa_params == "asm5") {
    smoothxg_poa_params = "-P 1,19,39,3,81,1"
  } else if (params.smoothxg_poa_params == "asm10") {
    smoothxg_poa_params = "-P 1,9,16,2,41,1"
  } else if (params.smoothxg_poa_params == "asm15") {
    smoothxg_poa_params = "-P 1,7,11,2,33,1"
  } else if (params.smoothxg_poa_params == "asm20") {
    smoothxg_poa_params = "-P 1,4,6,2,26,1"B
  } else {
    smoothxg_poa_params = "-P${params.smoothxg_poa_params}"
  }
}
def smoothxg_poa_params_display = smoothxg_poa_params.replaceAll(/,/, "_")
def smoothxg_temp_dir = params.smoothxg_temp_dir ? "-b${params.smoothxg_temp_dir}" : ""
def smoothxg_keep_intermediate_files = params.smoothxg_keep_intermediate_files ? "-K" : ""
def smoothxg_xpoa = "-S" 
if (params.smoothxg_run_abpoa != null) {
  smoothxg_xpoa = ""
}
def smoothxg_poa_mode = params.smoothxg_run_global_poa ? "-Z" : ""
// disabling consensus graph mode cause it is buggy
def smoothxg_consensus_spec = false

def wfmash_prefix = "wfmash"
def seqwish_prefix = ".seqwish"
def smoothxg_prefix = ".smoothxg"

def do_1d = true
def do_2d = true

if (params.no_viz) {
  do_1d = false
}

if (params.no_layout) {
  do_2d = false
}

def make_file_prefix = { f -> """\
${f.getName()}\
""" }

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


process wfmashMap {
  publishDir "${params.outdir}/wfmash_map", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta)
    path(fai)
    path(gzi)

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

process splitApproxMappingsInChunks {
  publishDir "${params.outdir}/wfmash_chunks", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(paf)
  output:
    path("${f}*.chunk_*.paf")
  """
  python3 /split_approx_mappings_in_chunks.py $paf ${params.wfmash_chunks}
  """
}

process wfmashAlign {
  publishDir "${params.outdir}/wfmash_align", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta), path(paf) 
    path(fai)
    path(gzi)

  output:
    path("${paf}.align.paf"), emit: paf

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
     -i $paf --invert-filtering \
     $fasta $fasta \
     >${paf}.align.paf
  """
}

process wfmash {
  publishDir "${params.outdir}/wfmash", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta)
    path(fai)
    path(gzi)

  output:
    tuple val(f), path("${f}.${wfmash_prefix}.paf")

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
     $fasta $fasta \
     >${f}.${wfmash_prefix}.paf
  """
}

process seqwish {
  publishDir "${params.outdir}/seqwish", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta)
    path(paf)

  output:
    tuple val(f), path("${f}${seqwish_prefix}.gfa")

  script:
    def input = paf.join(',')
    """
    seqwish \
      -t ${task.cpus} \
      -s $fasta \
      -p $input \
      -k ${params.seqwish_min_match_length} \
      -f ${params.seqwish_sparse_factor} \
      -g ${f}${seqwish_prefix}.gfa -P \
      -B ${params.seqwish_transclose_batch} \
      ${seqwish_temp_dir} \
      -P
    """
}

process smoothxg {
  publishDir "${params.outdir}/smoothxg", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(graph)

  output:
    path("${f}${smoothxg_prefix}.gfa"), emit: gfa_smooth
    path("${f}*.cons*.gfa"), optional: true, emit: consensus_smooth
    path("${f}.${smoothxg_prefix}.maf"), optional: true, emit: maf_smooth

  script:
    """
    maf_params=""
    if [[ ${params.smoothxg_write_maf} != false ]]; then
      maf_params="-m ${f}${smoothxg_prefix}.maf"
    fi
    smoothxg \
      -t ${task.cpus} \
      -T ${task.cpus} \
      -g ${graph} \
      -r ${n_haps} \
      ${smoothxg_temp_dir} \
      ${smoothxg_keep_intermediate_files} \
      -X 100 \
      -I ${smoothxg_block_id_min} \
      -R ${params.smoothxg_block_ratio_min} \
      -j ${params.smoothxg_max_path_jump} \
      -e ${params.smoothxg_max_edge_jump} \
      -l ${params.smoothxg_poa_length} \
      ${smoothxg_poa_params} \
      -O ${params.smoothxg_poa_padding} \
      -Y \$(echo "${params.smoothxg_pad_max_depth} * ${n_haps}" | bc) \
      -d 0 -D 0 \
      ${smoothxg_xpoa} \
      ${smoothxg_poa_mode} \
      \$maf_params \
      -V \
      -o ${f}${smoothxg_prefix}.gfa
    """
}

process gfaffix {
  publishDir "${params.outdir}/gfaffix", mode: "${params.publish_dir_mode}"

  input:
    path(graph)

  output:
    path("*.norm.og"), emit: og_norm
    path("${graph}.norm.gfa"), emit: gfa_norm
    path("${graph}.norm.affixes.tsv.gz"), emit: tsv_norm

  """
  gfaffix $graph -o ${graph}.norm.gfa | gzip > ${graph}.norm.affixes.tsv.gz 
  odgi build -g ${graph}.norm.gfa -o ${graph}.norm.gfa.og -P -t ${task.cpus} -O -o - \
  | odgi unchop -P -t ${task.cpus} -i - -o - \
  | odgi sort -i - -o ${graph}.norm.og -t ${task.cpus} -p Ygs 
  odgi view -i ${graph}.norm.og -g > ${graph}.norm.gfa
  """
}

process odgiBuild {
  publishDir "${params.outdir}/odgi_build", mode: "${params.publish_dir_mode}"

  input:
    path(graph)

  output:
    path("${graph}.og")

  """
  odgi build -g $graph -o ${graph}.og -P -t ${task.cpus}
  """
}

process odgiStats {
  publishDir "${params.outdir}/odgi_stats", mode: "${params.publish_dir_mode}"

  input:
    path(graph)

  output:
    path("${graph}.stats.yaml")

  """
  odgi stats -i "${graph}" -m > "${graph}.stats.yaml" 2>&1
  """
}

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

process odgiLayout {
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

process vg_deconstruct {
  publishDir "${params.outdir}/vg_deconstruct", mode: "${params.publish_dir_mode}"

  input:
  tuple path(graph), val(vcf_spec)

  output:
  path("${graph}.*.vcf"), emit: vg_deconstruct_vcf
  path("*.vcf.stats"), optional: true, emit: vg_deconstruct_bcftools_stats

  """
  ref=\$(echo "$vcf_spec" | cut -f 1 -d:)
  delim=\$(echo "$vcf_spec" | cut -f 2 -d:)
  pop_length=\$(echo "$vcf_spec" | cut -f 3 -d:)
  if [[ -z \$pop_length ]]; then
    pop_length=0
  fi
  vcf="${graph}".\$(echo \$ref | tr '/|' '_').vcf
  vg deconstruct -P \$ref -H \$delim -e -a -t "${task.cpus}" "${graph}" > \$vcf
  bcftools stats \$vcf > \$vcf.stats

  if [[ \$pop_length -gt 0 ]]; then
  vcf_decomposed=${graph}.final.\$(echo \$ref | tr '/|' '_').decomposed.vcf
  vcf_decomposed_tmp=\$vcf_decomposed.tmp.vcf
  bgzip -c -@ ${task.cpus} \$vcf > \$vcf.gz
  vcfbub -l 0 -a \$pop_length --input \$vcf.gz | vcfwave -I 1000 -t ${task.cpus} > \$vcf_decomposed_tmp

  #TODO: to remove when vcfwave will be bug-free
  # The TYPE info sometimes is wrong/missing
  # There are variants without the ALT allele
  bcftools annotate -x INFO/TYPE \$vcf_decomposed_tmp  | awk '\$5 != "."' > \$vcf_decomposed
  rm \$vcf_decomposed_tmp \$vcf.gz

  bcftools stats \$vcf_decomposed > \$vcf_decomposed.stats
fi
  """
}

// TODO ONCE OUR CUSTOM MULTIQC VERSION IS IN A MULTIQC RELEASE, WE CAN CHANGE THIS
process multiQC {
  publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

  input:
  path vg_deconstruct_bcftools_stats
  path odgi_stats
  path odgi_viz
  path odgi_draw
  path(multiqc_config)
  val(prefix)

  output:
  path "${prefix}_multiqc/*multiqc_report.html", emit: report
  path "${prefix}_multiqc/*_data"              , emit: data
  path "${prefix}_multiqc/*_plots"             , optional:true, emit: plots

  """
  multiqc -s . -c ${multiqc_config} --outdir ${prefix}_multiqc
  """
}

workflow PGGB {
  take:
  ch_fasta
  fai_path
  gzi_path

  main:

  fasta = ch_fasta.map { f -> tuple(make_file_prefix(f), f) }
  fasta_file_name = ch_fasta.map {it.getName()}
  
    if (!fai_path.exists() || !gzi_path.exists()) { // the assumption is that none of the files exist if only one does not exist
      samtoolsFaidx(fasta)
      fai = samtoolsFaidx.out.samtools_fai.collect()
      gzi = samtoolsFaidx.out.samtools_gzi.collect()
    } else {
      fai = channel.fromPath("${params.input}.fai").collect()
      gzi = channel.fromPath("${params.input}.gzi").collect()
    }
    if (params.wfmash_only) {
      // TODO Once we changed the way we changed the publish_dir_mode, we have to emit the .paf file as default, else not
      if (params.wfmash_chunks == 1) {
        wfmash(fasta, fai, gzi)
      } else {
        wfmashMap(fasta, fai, gzi)
        splitApproxMappingsInChunks(wfmashMap.out)
        // TODO update this once I understood it
        wfmashAlign(fasta.combine(splitApproxMappingsInChunks.out.flatten()), fai, gzi)
      }      
    } else {
      if (params.paf != null) {
        paf_ch = Channel.fromPath(params.paf)
        seqwish(fasta, paf_ch)
      } else {
        if (params.wfmash_chunks == 1) {
          wfmash(fasta, fai, gzi)
          seqwish(fasta, wfmash.out.collect{it[1]})
        } else {
          wfmashMap(fasta, fai, gzi)
          splitApproxMappingsInChunks(wfmashMap.out)
          wfmashAlign(fasta.combine(splitApproxMappingsInChunks.out.flatten()), fai, gzi)
          seqwish(fasta, wfmashAlign.out.collect())
        }
      }
      if (params.skip_smoothxg) {
        gfaffix(seqwish.out.collect{it[1]})  
        odgiBuild(seqwish.out.collect{it[1]})
      } else {
        smoothxg(seqwish.out)
        gfaffix(smoothxg.out.gfa_smooth)
        odgiBuild(seqwish.out.collect{it[1]}.mix(smoothxg.out.consensus_smooth.flatten()))
      }
      odgiStats(odgiBuild.out.mix(gfaffix.out.og_norm))

      odgiVizOut = Channel.empty()
      if (do_1d) {
          odgiVizOut = odgiViz(gfaffix.out.og_norm)
      }
      odgiDrawOut = Channel.empty()
      if (do_2d) {
        odgiLayout(gfaffix.out.og_norm)
        odgiDrawOut = odgiDraw(odgiLayout.out)
      }

      ch_vcf_spec = Channel.empty()
      vg_deconstruct = Channel.empty()
      if (params.vcf_spec != null) {
        ch_vcf_spec = Channel.from(params.vcf_spec).splitCsv().flatten()
        vg_deconstruct(gfaffix.out.gfa_norm.combine(ch_vcf_spec))
        // TODO add bcftools
        multiQC(
        vg_deconstruct.out.vg_deconstruct_bcftools_stats.collect().ifEmpty([]),
        odgiStats.out.collect().ifEmpty([]),
        odgiVizOut.collect().ifEmpty([]),
        odgiDrawOut.collect().ifEmpty([]),
        ch_multiqc_config
        )
      } else {
        multiQC(
          vg_deconstruct.collect().ifEmpty([]),
          odgiStats.out.collect().ifEmpty([]),
          odgiVizOut.collect().ifEmpty([]),
          odgiDrawOut.collect().ifEmpty([]),
          ch_multiqc_config,
          fasta_file_name
        )
      }
    }

    ch_empty = Channel.empty()
    // we currently don't want to emit anything

  emit:
    ch_empty

}