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

// We can't change global parameters inside this scope, so we build the ones we need locally
def alignment_merge_cmd = params.alignment_merge_segments ? "-M" : params.alignment_merge_cmd
def alignment_exclude_cmd = params.alignment_exclude_delim ? "-Y${params.alignment_exclude_delim}" : params.alignment_exclude_cmd
def alignment_split_cmd = params.alignment_no_splits ? "-N" : params.alignment_split_cmd
def aligner = params.wfmash ? "W" : "E"
def edyeet_align_pct_id_display = params.wfmash ? "" : "a${params.edyeet_align_pct_id}-"
def smoothxg_poa_params_display = params.smoothxg_poa_params.replaceAll(/,/, "_")
def file_name_prefix_display = ""
def alignment_prefix = ""
def seqwish_prefix = ""
def smoothxg_prefix = ""
// default case
if (!params.file_name_prefix) {
  file_name_prefix_display = ".pggb"
  alignment_prefix = "-${aligner}"
  seqwish_prefix = ".seqwish"
  smoothxg_prefix = ".smoothxg"
} else if (params.file_name_prefix == "pggb") {
  // fancy naming scheme
  file_name_prefix_display = ".pggb"
  alignment_prefix = """-\
${aligner}-\
s${params.alignment_segment_length}-\
l${params.alignment_block_length}-\
p${params.alignment_map_pct_id}-\
n${params.alignment_n_secondary}-\
${edyeet_align_pct_id_display}\
K${params.alignment_mash_kmer}\
${alignment_merge_cmd}\
${alignment_split_cmd}\
${alignment_exclude_cmd}\
"""
  seqwish_prefix = """\
.seqwish-\
k${params.seqwish_min_match_length}-\
B${params.seqwish_transclose_batch}\
"""
  smoothxg_prefix = """${seqwish_prefix}\
.smoothxg-\
w${params.smoothxg_max_block_weight}-\
j${params.smoothxg_max_path_jump}-\
e${params.smoothxg_max_edge_jump}-\
I${params.smoothxg_block_id_min}-\
p${smoothxg_poa_params_display}-M-J0.7-K-G150\
"""
} else {
  // take the given prefix
  file_name_prefix_display= "${params.file_name_prefix}.pggb"
  aligment_prefix = "-${aligner}"
  seqwish_prefix = ".seqwish"
  smoothxg_prefix = ".smoothxg"
}

// TODO update, make_file_prefix
def make_file_prefix = { f -> """\
${f.getName()}${file_name_prefix_display}\
""" }

def make_file_prefix_given = { f -> """\
${file_name_prefix_display}\
""" }

if (!params.file_name_prefix || params.file_name_prefix == "pggb") {
  fasta = channel.fromPath("${params.input}").map { f -> tuple(make_file_prefix(f), f) }
} else {
  fasta = channel.fromPath("${params.input}").map { f -> tuple(make_file_prefix_given(f), f) }
}

process edyeet {
  input:
  tuple val(f), path(fasta)

  output:
  tuple val(f), path(fasta), path("${f}${alignment_prefix}.paf")

  """
  edyeet ${alignment_exclude_cmd} \
     -s ${params.alignment_segment_length} \
     -l ${params.alignment_block_length} \
     ${alignment_merge_cmd} \
     ${alignment_split_cmd} \
     -p ${params.alignment_map_pct_id} \
     -n ${params.alignment_n_secondary} \
     -a ${params.edyeet_align_pct_id} \
     -k ${params.alignment_mash_kmer} \
     -t ${task.cpus} \
     $fasta $fasta \
     >${f}${alignment_prefix}.paf 
  """
}

process wfmash {
  input:
  tuple val(f), path(fasta)

  output:
  tuple val(f), path(fasta), path("${f}${alignment_prefix}.paf")

  """
  wfmash ${alignment_exclude_cmd} \
     -s ${params.alignment_segment_length} \
     -l ${params.alignment_block_length} \
     ${alignment_merge_cmd} \
     ${alignment_split_cmd} \
     -p ${params.alignment_map_pct_id} \
     -n ${params.alignment_n_secondary} \
     -k ${params.alignment_mash_kmer} \
     -t ${task.cpus} \
     $fasta $fasta \
     >${f}${alignment_prefix}.paf 
  """
}

process seqwish {
  publishDir "${params.outdir}/seqwish", mode: "${params.publish_dir_mode}"

  input:
  tuple val(f), path(fasta), path(alignment)

  output:
    tuple val(f), path("${f}${seqwish_prefix}.gfa")

  script:
    """
    seqwish \
      -t ${task.cpus} \
      -s $fasta \
      -p $alignment \
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
    path("${f}*.consensus*.gfa"), emit: consensus_smooth
    path("${f}${smoothxg_prefix}.maf"), emit: maf_smooth

  script:
    """
    smoothxg \
      -t ${task.cpus} \
      -g $graph \
      -w ${params.smoothxg_max_block_weight} \
      -M \
      -J 0.7 \
      -K \
      -G 150 \
      -I ${params.smoothxg_block_id_min} \
      -R ${params.smoothxg_ratio_contain} \
      -j ${params.smoothxg_max_path_jump} \
      -e ${params.smoothxg_max_edge_jump} \
      -l ${params.smoothxg_max_poa_length} \
      -p ${params.smoothxg_poa_params} \
      -m ${f}${smoothxg_prefix}.maf \
      -C ${f}${smoothxg_prefix}.consensus,${params.smoothxg_consensus_spec} \
      -o ${f}${smoothxg_prefix}.gfa
    """  
}

process odgiBuild {
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
  path("${graph}.stats")

  """
  odgi stats -i "${graph}" -S -s -d -l > "${graph}.stats" 2>&1
  """
}

process odgiViz {
  publishDir "${params.outdir}/odgi_viz", mode: "${params.publish_dir_mode}"

  input:
  path(graph)

  output:
  path("${graph}.viz_mqc.png")

  """
  odgi viz \
    -i $graph \
    -o ${graph}.viz_mqc.png \
    -x 1500 -y 500 -P 5
  """
}

process odgiChop {
  input:
  path(graph)

  output:
  path("${graph}.chop.og")

  """
  odgi chop -i $graph -c 100 -o ${graph}.chop.og
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
  path("${layoutGraph}.draw_mqc.png")

  """
  odgi draw \
    -i $graph \
    -c $layoutGraph \
    -p ${layoutGraph}.draw_mqc.png \
    -C \
    -w 50 \
    -H 1000 -t ${task.cpus}
  """
}

process multiQC {
  publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

  input:
  path odgi_stats
  path odgi_viz
  path odgi_draw

  output:
  path "*multiqc_report.html", emit: report
  path "*_data"              , emit: data
  path "*_plots"             , optional:true, emit: plots

  """
  multiqc -s .
  """
}

workflow {
  main:
    if (params.wfmash == false) {
      edyeet(fasta)
      seqwish(edyeet.out)
    } else {
      wfmash(fasta)
      seqwish(wfmash.out)
    }
    smoothxg(seqwish.out)
    if (params.do_stats) { 
      odgiBuild(seqwish.out.collect{it[1]}.mix(smoothxg.out.gfa_smooth, smoothxg.out.consensus_smooth.flatten())) 
      odgiStats(odgiBuild.out)  
    }
    else {
      odgiBuild(smoothxg.out.gfa_smooth)
    }
    odgiVizOut = Channel.empty()
    if (params.do_viz) {
      odgiVizOut = odgiViz(odgiBuild.out.filter(~".*smooth.*"))
    }
    odgiDrawOut = Channel.empty()
    if (params.do_layout) {
      odgiChop(odgiBuild.out.filter(~".*smooth.*"))
      odgiLayout(odgiChop.out)
      odgiDrawOut = odgiDraw(odgiLayout.out)
    }

    multiQC(
      odgiStats.out.collect().ifEmpty([]),
      odgiVizOut.collect().ifEmpty([]),
      odgiDrawOut.collect().ifEmpty([])
    )
}

// /*
//  * Stage config files
//  */
// ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
// ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
// ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

// /*
//  * Check parameters
//  */
// Checks.aws_batch(workflow, params)     // Check AWS batch settings
// Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

// /*
//  * Print parameter summary
//  */
// // Has the run name been specified by the user?
// // this has the bonus effect of catching both -name and --name
// run_name = params.name
// if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
//     run_name = workflow.runName
// }
// summary = Schema.params_summary(workflow, params, run_name)
// log.info Headers.nf_core(workflow, params.monochrome_logs)
// log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
// log.info "-\033[2m----------------------------------------------------\033[0m-"

// workflow_summary = Schema.params_mqc_summary(summary)
// ch_workflow_summary = Channel.value(workflow_summary)

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

// /*
//  * Run the workflow
//  */
// workflow {

//     CHECK_SAMPLESHEET(ch_input)
//         .splitCsv(header:true, sep:',')
//         .map { check_samplesheet_paths(it) }
//         .set { ch_raw_reads }

//     FASTQC(ch_raw_reads)

//     OUTPUT_DOCUMENTATION(
//         ch_output_docs,
//         ch_output_docs_images)

//     GET_SOFTWARE_VERSIONS()

//     MULTIQC(
//         ch_multiqc_config,
//         ch_multiqc_custom_config.collect().ifEmpty([]),
//         FASTQC.out.collect(),
//         GET_SOFTWARE_VERSIONS.out.yml.collect(),
//         ch_workflow_summary)
// }

// /*
//  * Send completion email
//  */
// workflow.onComplete {
//     def multiqc_report = []
//     Completion.email(workflow, params, summary, run_name, baseDir, multiqc_report, log)
//     Completion.summary(workflow, params, log)
// }

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
    Alignment options:
      --edyeet_align_pct_id [n]       percent identity in the edyeet edlib alignment step [default: 90]
      --alignment_map_pct_id [n]      percent identity in the wfmash or edyeet mashmap step [default: 90]
      --alignment_n_secondary [n]     number of secondary mappings to retain in 'map' filter mode [default: 10]
      --alignment_segment_length [n]  segment length for mapping [default: 10000]
      --alignment_block_length [n]    minimum block length filter for mapping [default: 3*alignment_segment_length]
      --alignment_mash_kmer [n]       kmer size for mashmap [default: 16]
      --alignment_merge_segments      merge successive mappings [default: OFF]
      --alignment_no_splits           disable splitting of input sequences during mapping [default: OFF]
      --alignment_exclude--delim [c]  skip mappings between sequences with the same name prefix before
                                      the given delimiter character [default: all-vs-all and !self]
    Seqwish options:
      --seqwish_min_match_length [n]  ignore exact matches below this length [default: 19]
      --seqwish_transclose_batch [n]  number of bp to use for transitive closure batch [default: 1000000]

    Smoothxg options:
      --smoothxg_max_block_weight [n] maximum seed sequence in block [default: 10000]
      --smoothxg_max_path_jump [n]    maximum path jump to include in block [default: 5000]
      --smoothxg_max_edge_jump [n]    maximum edge jump before breaking [default: 5000]
      --smoothxg_max_popa_length [n]  maximum sequence length to put into POA [default: 10000]
      --smoothxg_consensus_spec [str] consensus graph specification: write the consensus graph to
                                      BASENAME.cons_[spec].gfa; where each spec contains at least a min_len parameter
                                      (which defines the length of divergences from consensus paths to preserve in the
                                      output), optionally a file containing reference paths to preserve in the output,
                                      a flag (y/n) indicating whether we should also use the POA consensus paths, a
                                      minimum coverage of consensus paths to retain (min_cov), and a maximum allele
                                      length (max_len, defaults to 1e6); implies -a; example:
                                      cons,100,1000:refs1.txt:n,1000:refs2.txt:y:2.3:1000000,10000
                                      [default: 10,100,1000,10000]
      --smoothxg_block_id_min [n]     split blocks into groups connected by this identity threshold [default: OFF]
      --smoothxg_ratio_contain [n]    minimum short length / long length ratio to compare sequences for the containment
                                      metric in the clustering [default: 0]
      --smoothxg_poa_params [str]     score parameters for POA in the form of match,mismatch,gap1,ext1,gap2,ext2
                                      [default: 1,4,6,2,26,1]                                         

    Visualization options:
      --do_viz                        Generate 1D visualisations of the built graphs [default: OFF]
      --do_layout                     Generate 2D visualisations of the built graphs [default: OFF]

    Other options:
      --outdir [file]                 The output directory where the results will be saved [default: ./results]
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
      --file_name_prefix [str]        Prefix for the output file names. If 'pggb', the file names will be very verbose and contain all parameters for each process. [default: --input]

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
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
