#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/pangenome
========================================================================================
 nf-core/pangenome Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/pangenome
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

if (params.input == null) {
        log.info"""

    Mandatory argument --input missing! For more details run with --help.

    """.stripIndent()  

    exit 1
}

if (params.n_mappings == null) {
        log.info"""

    Mandatory argument --wfmash_n_mappings missing! For more details run with --help.

    """.stripIndent()  

    exit 1
}

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()

// We can't change global parameters inside this scope, so we build the ones we need locally
def n_haps = 0
if (!params.smoothxg_num_haps) {
  n_haps = params.n_mappings
}

def wfmash_merge_cmd = params.wfmash_merge_segments ? "-M" : ""
def wfmash_exclude_cmd = params.wfmash_exclude_delim ? "-Y${params.wfmash_exclude_delim}" : "-X"
def wfmash_split_cmd = params.wfmash_no_splits ? "-N" : ""
def wfmash_block_length_cmd = params.wfmash_block_length ? "-l${params.wfmash_block_length}" : ""
def wfmash_mash_kmer_cmd = params.wfmash_mash_kmer ? "-k${params.wfmash_mash_kmer}" : ""
def wfmash_kmer_thres_cmd = params.wfmash_mash_kmer_thres ? "-H${params.wfmash_kmer_thres}" : ""
def wfmash_n_mappings_minus_1 = params.n_mappings - 1
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

def smoothxg_block_id_min = params.wfmash_map_pct_id / 100.0
// TODO: CHANGE TO LARGE P ONCE WE ARE THERE
def smoothxg_poa_params_cmd = ""
if (params.smoothxg_poa_params == null) {
  smoothxg_poa_params = "-p 1,19,39,3,81,1"
} else {
  if (params.smoothxg_poa_params == "asm5") {
    smoothxg_poa_params = "-p 1,19,39,3,81,1"
  } else if (params.smoothxg_poa_params == "asm10") {
    smoothxg_poa_params = "-p 1,9,16,2,41,1"
  } else if (params.smoothxg_poa_params == "asm15") {
    smoothxg_poa_params = "-p 1,7,11,2,33,1"
  } else if (params.smoothxg_poa_params == "asm20") {
    smoothxg_poa_params = "-p 1,4,6,2,26,1"
  } else {
    smoothxg_poa_params = "-p${params.smoothxg_poa_params}"
  }
}
def smoothxg_poa_params_display = smoothxg_poa_params.replaceAll(/,/, "_")
def smoothxg_temp_dir = params.smoothxg_temp_dir ? "-b${params.smoothxg_temp_dir}" : ""
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

fasta = channel.fromPath("${params.input}").map { f -> tuple(make_file_prefix(f), f) }

process wfmashMap {
  publishDir "${params.outdir}/wfmash_map", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta)

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

  output:
    path("${paf}.align.paf")

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
     -t ${task.cpus} \
     -i $paf \
     $fasta $fasta \
     >${paf}.align.paf
  """
}

process wfmash {
  publishDir "${params.outdir}/wfmash", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta)

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
     -t ${task.cpus} \
     $fasta $fasta \
     >${f}.${wfmash_prefix}.paf
  """
}

/// TODO -f ${params.seqwish_sparse_factor} \
process seqwish {
  publishDir "${params.outdir}/seqwish", mode: "${params.publish_dir_mode}"

  input:
    tuple val(f), path(fasta)
    path(pafs)

  output:
    tuple val(f), path("${f}${seqwish_prefix}.gfa")

  script:
    """
    if [[ \$(ls *.paf | wc -l) == 1 ]]; then
      input=$pafs
    else 
      input=\$(ls *.paf | tr '\\\n' ',')
      input=\${input::-1}
    fi  
    seqwish \
      -t ${task.cpus} \
      -s $fasta \
      -p \$input \
      -k ${params.seqwish_min_match_length} \
      -g ${f}${seqwish_prefix}.gfa -P \
      -B ${params.seqwish_transclose_batch} \
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
    smooth_iterations=\$(echo ${params.smoothxg_poa_length} | tr ',' '\\\n' | wc -l)
    echo \$smooth_iterations > smooth_iterations
    maf_params=""
    if [[ ${params.smoothxg_write_maf} != false ]]; then
      maf_params="-m ${f}.${smoothxg_prefix}.maf"
    fi
    for i in \$(seq 1 \$smooth_iterations);
    do
      input_gfa=${graph}
      if [[ \$i != 1 ]]; then
        input_gfa=smooth.\$(echo \$i - 1 | bc).gfa 
      fi
      if [[ \$i != \$smooth_iterations ]]; then
        poa_length=\$(echo ${params.smoothxg_poa_length} | cut -f \$i -d,)
        smoothxg \
          -t ${task.cpus} \
          -T ${task.cpus} \
          -g \$input_gfa \
          -w \$(echo "\$poa_length * ${n_haps}" | bc) \
          ${smoothxg_temp_dir} \
          -K \
          -X 100 \
          -I ${smoothxg_block_id_min} \
          -R ${params.smoothxg_block_ratio_min} \
          -j ${params.smoothxg_max_path_jump} \
          -e ${params.smoothxg_max_edge_jump} \
          -l \$poa_length \
          ${smoothxg_poa_params} \
          -O ${params.smoothxg_poa_padding} \
          -Y \$(echo "${params.smoothxg_pad_max_depth} * ${n_haps}" | bc) \
          -d 0 -D 0 \
          -V \
          -o smooth.\$i.gfa
      else
        poa_length=\$(echo ${params.smoothxg_poa_length} | cut -f \$i -d,)
        consensus_params="-V"
        if [[ ${params.smoothxg_consensus_spec} != false ]]; then
          consensus_params="-C ${f}.cons,${params.smoothxg_consensus_spec}"  
        fi  
        smoothxg \
          -t ${task.cpus} \
          -T ${task.cpus} \
          -g \$input_gfa \
          -w \$(echo "\$poa_length * ${n_haps}" | bc) \
          ${smoothxg_temp_dir} \
          -K \
          -X 100 \
          -I ${smoothxg_block_id_min} \
          -R ${params.smoothxg_block_ratio_min} \
          -j ${params.smoothxg_max_path_jump} \
          -e ${params.smoothxg_max_edge_jump} \
          -l \$poa_length \
          ${smoothxg_poa_params} \
          -O ${params.smoothxg_poa_padding} \
          -Y \$(echo "${params.smoothxg_pad_max_depth} * ${n_haps}" | bc) \
          -d 0 -D 0 \
          \$maf_params \
          -Q ${params.smoothxg_consensus_prefix} \
          \$consensus_params \
          -o ${f}${smoothxg_prefix}.gfa
      fi
    done  
    """
}

process gfaffix {
  publishDir "${params.outdir}/gfaffix", mode: "${params.publish_dir_mode}"

  input:
    path(graph)

  output:
    path("${graph}.norm.og"), emit: og_norm
    path("${graph}.norm.gfa"), emit: gfa_norm
    path("${graph}.norm.affixes.tsv.gz"), emit: tsv_norm

  """
  gfaffix $graph -o ${graph}.norm.gfa | gzip > ${graph}.norm.affixes.tsv.gz 
  odgi build -g ${graph}.norm.gfa -o ${graph}.norm.gfa.og -P -t ${task.cpus} -O -o - \
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

// TODO we can parallelize this for each reference given in ${params.vcf_spec}
process vg_deconstruct {
  publishDir "${params.outdir}/vg_deconstruct", mode: "${params.publish_dir_mode}"

  input:
  tuple path(graph), val(vcf_spec)

  output:
  path("${graph}.*.vcf")

  """
  ref=\$(echo "$vcf_spec" | cut -f 1 -d:)
  delim=\$(echo "$vcf_spec" | cut -f 2 -d:)
  vcf="${graph}".\$(echo \$ref | tr '/|' '_').vcf
  vg deconstruct -P \$ref -H \$delim -e -a -t "${task.cpus}" "${graph}" > \$vcf
  """
}

// TODO ONCE OUR CUSTOM MULTIQC VERSION IS IN A MULTIQC RELEASE, WE CAN CHANGE THIS
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

workflow {
  main:

    if (params.wfmash_only) {
      // TODO Once we changed the way we changed the publish_dir_mode, we have to emit the .paf file as default, else not
      if (params.wfmash_chunks == 1) {
        wfmash(fasta)
      } else {
        wfmashMap(fasta)
        splitApproxMappingsInChunks(wfmashMap.out)
        wfmashAlign(fasta.combine(splitApproxMappingsInChunks.out.flatten()))
      }      
    } else {
      if (params.paf != false) {
        paf_ch = Channel.fromPath(params.paf)
        seqwish(fasta, paf_ch)
      } else {
        if (params.wfmash_chunks == 1) {
          wfmash(fasta)
          seqwish(fasta, wfmash.out.collect{it[1]})
        } else {
          wfmashMap(fasta)
          splitApproxMappingsInChunks(wfmashMap.out)
          wfmashAlign(fasta.combine(splitApproxMappingsInChunks.out.flatten()))
          seqwish(fasta, wfmashAlign.out.collect())
        }
      }
      smoothxg(seqwish.out)
      gfaffix(smoothxg.out.gfa_smooth)

      odgiBuild(seqwish.out.collect{it[1]}.mix(smoothxg.out.consensus_smooth.flatten()))
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

      ch_vg_deconstruct = Channel.empty()
      ch_vcf_spec = Channel.empty()
      if (params.vcf_spec != false) {
        ch_vcf_spec = Channel.from(params.vcf_spec).splitCsv().flatten()
        ch_vg_deconstruct = vg_deconstruct(gfaffix.out.gfa_norm.combine(ch_vcf_spec))
        // TODO add bcftools
      }

      multiQC(
        odgiStats.out.collect().ifEmpty([]),
        odgiVizOut.collect().ifEmpty([]),
        odgiDrawOut.collect().ifEmpty([]),
        ch_multiqc_config
      )
    }
}

// /*
//  * Include local pipeline modules
//  */
// include { OUTPUT_DOCUMENTATION } from './modules/local/output_documentation' params(params)
// include { GET_SOFTWARE_VERSIONS } from './modules/local/get_software_versions' params(params)
// include { CHECK_SAMPLESHEET; check_samplesheet_paths } from './modules/local/check_samplesheet' params(params)

// /*
//  * Include nf-core modules
//  */
// include { FASTQC } from './modules/nf-core/fastqc' params(params)
// include { MULTIQC } from './modules/nf-core/multiqc' params(params)

// TODO REMOVE THIS ONCE WE TOOK CARE OF nextflow_schema.json
def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/pangenome --input 'data/input.fa.gz' -profile docker

    Mandatory arguments:
      --input [file]                  Path to input FASTA (must be surrounded with quotes)
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more
    Wfmash options:
      --wfmash_map_pct_id [n]         percent identity in the wfmash mashmap step [default: 90]
      --wfmash_n_mappings [n]         number of secondary mappings to retain in 'map' filter mode [default: 10]
      --wfmash_segment_length [n]     segment length for mapping [default: 3000]
      --wfmash_block_length [n]       minimum block length filter for mapping [default: 3 * wfmash_segment_length]
      --wfmash_mash_kmer [n]          kmer size for mashmap [default: 16]
      --wfmash_merge_segments         merge successive mappings [default: OFF]
      --wfmash_no_splits              disable splitting of input sequences during mapping [default: OFF]
      --wfmash_exclude--delim [c]     skip mappings between sequences with the same name prefix before
                                      the given delimiter character [default: all-vs-all and !self]
      --wfmash_chunks                 The number of files to generate from the approximate wfmash mappings to scale across a whole cluster. It is recommended to set this to the number of available nodes. If only one machine is available, leave it at 1. [default: 1]
      --wfmash_only                   If this parameter is set, only the wfmash alignment step of the pipeline is executed. This option is offered for users who want to use wfmash on a cluster. [default: OFF]

    Seqwish options:
      --seqwish_min_match_length [n]  ignore exact matches below this length [default: 47]
      --seqwish_transclose_batch [n]  number of bp to use for transitive closure batch [default: 10000000]

    Smoothxg options:
      --smoothxg_num_haps [n]         number of haplotypes in the given FASTA [default: wfmash_n_mappings]
      --smoothxg_max_path_jump [n]    maximum path jump to include in block [default: 0]
      --smoothxg_max_edge_jump [n]    maximum edge jump before breaking [default: 0]
      --smoothxg_max_poa_length [n]   maximum sequence length to put into POA, can be a comma-separated list; 
                                      for each element smoothxg will be executed once [default: 4001,4507]
      --smoothxg_consensus_spec [str] consensus graph specification: write the consensus graph to
                                      BASENAME.cons_[spec].gfa; where each spec contains at least a min_len parameter
                                      (which defines the length of divergences from consensus paths to preserve in the
                                      output), optionally a file containing reference paths to preserve in the output,
                                      a flag (y/n) indicating whether we should also use the POA consensus paths, a
                                      minimum coverage of consensus paths to retain (min_cov), and a maximum allele
                                      length (max_len, defaults to 1e6); implies -a; example:
                                      cons,100,1000:refs1.txt:n,1000:refs2.txt:y:2.3:1000000,10000
                                      [default: OFF]
      --smoothxg_consensus_prefix [n] use this prefix for consensus path names [default: Consensus_]
      --smoothxg_block_ratio_min [n]  minimum small / large length ratio to cluster in a block [default: 0.0]
      --smoothxg_block_id_min [n]     split blocks into groups connected by this identity threshold [default: 0.95]
      --smoothxg_pad_max_depth [n]    path depth at which we don't pad the POA problem [default: 100]
      --smoothxg_poa_padding [n]      pad each end of each sequence in POA with N*(longest_poas_seq) bp [default: 0.03]
      --smoothxg_poa_params [str]     score parameters for POA in the form of match,mismatch,gap1,ext1,gap2,ext2
                                      [default: 1,19,39,3,81,1]
      --smoothxg_write_maf [n]        write MAF output representing merged POA blocks [default: OFF]

    Visualization options:
      --viz                           Generate 1D and 2D visualisations of the built graphs [default: OFF]

    VCF options:
      --vcf_spec                      specify a set of VCFs to produce with SPEC = REF:DELIM[,REF:DELIM]*
                                      the paths matching ^REF are used as a reference, while the sample haplotypes
                                      are derived from path names, e.g. when DELIM=# and with '-V chm13:#',
                                      a path named HG002#1#ctg would be assigned to sample HG002 phase 1 [default: OFF]

    Other options:
      --outdir [file]                 The output directory where the results will be saved [default: ./results]
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
// TODO INVOKE THIS AGAIN ONCE IT IS CLEAR HOW TO ADD A NAME TO THE RUN
// TODO ERROR: You used a core Nextflow option with two hyphens: '--name'. Please resubmit with '-name'
/*
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}
*/

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+

// TODO ACTIVATE AFTER nextflow_schema.json has bin updated

def json_schema = "$projectDir/nextflow_schema.json"
/*
if (params.help) {
    def command = "nextflow run nf-core/pangenome --input "input.fasta.gz" -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}
*/

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+

// TODO ACTIVATE AFTER nextflow_schema.json has bin updated
/*
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}
*/

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// TODO UNTICK THIS WHEN WE ACTUALLY USE IT
// Check if genome exists in the config file
/*
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
}
*/

// TODO nf-core: Add any reference files that are needed
// Configurable reference genomes
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the channel below in a process, define the following:
//   input:
//   file fasta from ch_fasta
//

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
// ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)


////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']            = params.input
summary['Data Type']        = 'FASTA'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}

// TODO DEACTIVATE AFTER nextflow_schema.json has bin updated
def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/pangenome v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
// TODO DEACTIVATE AFTER nextflow_schema.json has bin updated

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-pangenome-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/pangenome Workflow Summary'
    section_href: 'https://github.com/nf-core/pangenome'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 1 - FastQC
 */
 /// THIS IS NOT NEEDED IN THIS PIPELINE
 /*
process fastqc {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf('.zip') > 0 ? "zips/$filename" : "$filename"
        }

    input:
    set val(name), file(reads) from ch_read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc_results

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}
*/

/*
 * STEP 2 - MultiQC
 */
 /// we currently have to use our own version of MultiQC that is not release, yet
 /*
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    // TODO nf-core: Add in log files from your new processes for MultiQC to find!
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}
*/

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
      file output_docs from ch_output_docs
      file images from ch_output_docs_images

    output:
      file 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/pangenome] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/pangenome] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/pangenome] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/pangenome] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/pangenome] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/pangenome] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/pangenome]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/pangenome]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error '====================================================\n' +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            '============================================================'
                }
            }
        }
    }
}
