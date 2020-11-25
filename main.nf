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

params.input="${baseDir}/**.fa.gz"
params.map_pct_id=90//false
params.align_pct_id=90// false
params.n_secondary=10//false
params.segment_length=1000//false
params.mash_kmer=16
params.min_match_length=8
params.max_block_weight=10000
params.max_path_jump=5000
params.min_subpath=0
params.max_edge_jump=5000
params.max_poa_length=10000
params.do_viz=false
params.do_layout=false
params.do_stats=false

def makeBaseName = { f -> """\
${f.getSimpleName()}.pggb-\
s${params.segment_length}-\
p${params.map_pct_id}-\
n${params.n_secondary}-\
a${params.align_pct_id}-\
K${params.mash_kmer}-\
k${params.min_match_length}-\
w${params.max_block_weight}-\
j${params.max_path_jump}-\
W${params.min_subpath}-\
e${params.max_edge_jump}\
""" }

fasta = channel.fromPath("${params.input}").map { f -> tuple(makeBaseName(f), f) }

process edyeet {
  input:
  tuple val(f), path(fasta)

  output:
  tuple val(f), path(fasta), path("${f}.paf")

  """
  edyeet -X \
     -s ${params.segment_length} \
     -p ${params.map_pct_id} \
     -n ${params.n_secondary} \
     -a ${params.align_pct_id} \
     -k ${params.mash_kmer} \
     -t ${task.cpus} \
     $fasta $fasta \
     >${f}.paf 
  """
}


process seqwish {
  input:
  tuple val(f), path(fasta), path(alignment)

  output:
  tuple val(f), path("${f}.seqwish.gfa")

  """
  seqwish \
    -t ${task.cpus} \
    -s $fasta \
    -p $alignment \
    -k ${params.min_match_length} \
    -g ${f}.seqwish.gfa -P
  """
}

process smoothxg {
  input:
    tuple val(f), path(graph)

  output:
    // val(f), emit: failure
    path("${f}.smooth.gfa"), emit: gfa_smooth // TODO, file("${f}*.consensus*.gfa")

  script:
    """
    smoothxg \
      -t ${task.cpus} \
      -g $graph \
      -w ${params.max_block_weight} \
      -j ${params.max_path_jump} \
      -e ${params.max_edge_jump} \
      -l ${params.max_poa_length} \
      -o ${f}.smooth.gfa \
      -m ${f}.smooth.maf \
      -s ${f}.consensus \
      -a \
      -C 5000
    """  
}

process odgiBuild {
  input:
  // tuple val(f), file(graph), file(consensus_graphs)
  path(graph)

  output:
  path("${graph}.og")

  """
  odgi build -g $graph -o ${graph}.og
  """
}

process odgiStats {
  publishDir "${params.outdir}/odgiStats", mode: 'copy'

  input: 
  path(graph)

  output:
  path("${graph}.stats")

  """
  odgi stats -i $graph -S -s -d -l > "${graph}.stats"
  """
}

process odgiViz {
  input:
  path(graph)

  output:
  path("${graph}.viz.png")

  """
  odgi viz \
    -i $graph \
    -o ${graph}.viz.png \
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
  input:
  tuple path(graph), path(layoutGraph)

  output:
  path("${layoutGraph}.png")

  """
  odgi draw \
    -i $graph \
    -c $layoutGraph \
    -p ${layoutGraph}.png \
    -H 1000 -t ${task.cpus}
  """
}

workflow {
    edyeet(fasta)
    seqwish(edyeet.out)
    smoothxg(seqwish.out)
    odgiBuild(smoothxg.out.gfa_smooth)
    odgiStats(odgiBuild.out)
    odgiViz(odgiBuild.out)
    odgiChop(odgiBuild.out)
    odgiLayout(odgiChop.out)
    odgiDraw(odgiLayout.out)
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

    nextflow run nf-core/pangenome --input 'data/*.fa.gz' -profile docker

    Mandatory arguments:
      --input [file]                  Path to input data (must be surrounded with quotes)
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

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
