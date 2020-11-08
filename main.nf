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

nextflow.preview.dsl = 2

/*
 * Print help message if required
 */
if (params.help) {
    // TODO nf-core: Update typical command used to run pipeline
    def command = "nextflow run nf-core/pangenome --input samplesheet.csv -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$baseDir/nextflow_schema.json", command)
    exit 0
}

/*
* TODO integrate this later....
*/

/*
params.input_fasta="${baseDir}/**.fa.gz"
params.map_pct_id=false
params.align_pct_id=false
params.n_secondary=false
params.segment_length=false
params.mash_kmer=16
params.min_match_length=8
params.max_block_weight=10000
params.max_path_jump=5000
params.min_subpath=0
params.max_edge_jump=5000
params.max_poa_length=10000
params.do_viz=false
params.do_layout=false

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

fasta = Channel.fromPath("${params.input_fasta}").map { f -> tuple(makeBaseName(f), f) }

process edyeet {
  input:
  tuple val(f), file(fasta)

  output:
  tuple val(f), file(fasta), file("${f}.paf")

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
  tuple val(f), file(fasta), file(alignment)

  output:
  tuple val(f), file("${f}.seqwish.gfa")

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
  tuple val(f), file(graph)

  output:
  tuple val(f), file("${f}.smooth.gfa")

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
    -C 5000 \
  """
}

process odgiBuild {
  input:
  tuple val(f), file(graph)

  output:
  tuple val(f), file("${f}.smooth.og")

  """
  odgi build -g $graph -o ${f}.smooth.og
  """
}

process odgiViz {
  input:
  tuple val(f), file(graph)

  output:
  tuple val(f), file("${f}.smooth.og.viz.png")

  """
  odgi viz \
    -i $graph \
    -o ${f}.smooth.og.viz.png \
    -x 1500 -y 500 -P 5
  """
}

process odgiChop {
  input:
  tuple val(f), file(graph)

  output:
  tuple val(f), file("${f}.smooth.chop.og")

  """
  odgi chop -i $graph -c 100 -o ${f}.smooth.chop.og
  """
}

process odgiLayout {
  input:
  tuple val(f), file(graph)

  output:
  tuple val(f), file(graph), file("${f}.smooth.chop.og.lay")

  """
  odgi layout \
    -i $graph \
    -o ${f}.smooth.chop.og.lay \
    -t ${task.cpus} -P
  """
}

process odgiDraw {
  input:
  tuple val(f), file(graph), file(layoutGraph)

  output:
  tuple val(f), file("${f}.smooth.chop.og.lay.png")

  """
  odgi draw \
    -i $graph \
    -c $layoutGraph \
    -p ${f}.smooth.chop.og.lay.png \
    -H 1000 -t ${task.cpus}
  """
}


workflow {
    edyeet(fasta)
    seqwish(edyeet.out)
    smoothxg(seqwish.out)
    odgiBuild(smoothxg.out)
    odgiViz(odgiBuild.out)
    odgiChop(odgiBuild.out)
    odgiLayout(odgiChop.out)
    odgiDraw(odgiLayout.out)
}
*/

/*
 * Stage config files
 */
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

/*
 * Validate parameters
 */
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }

/*
 * Reference genomes
 */
// TODO nf-core: Add any reference files that are needed
// NOTE - FOR SIMPLICITY THIS IS NOT USED IN THIS PIPELINE
// EXAMPLE ONLY TO DEMONSTRATE USAGE OF AWS IGENOMES
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}
params.fasta = params.genomes[params.genome]?.fasta
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) }

/*
 * Check parameters
 */
Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

/*
 * Print parameter summary
 */
// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    run_name = workflow.runName
}
summary = Schema.params_summary(workflow, params, run_name)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

workflow_summary = Schema.params_mqc_summary(summary)
ch_workflow_summary = Channel.value(workflow_summary)

/*
 * Include local pipeline modules
 */
include { OUTPUT_DOCUMENTATION } from './modules/local/output_documentation' params(params)
include { GET_SOFTWARE_VERSIONS } from './modules/local/get_software_versions' params(params)
include { CHECK_SAMPLESHEET; check_samplesheet_paths } from './modules/local/check_samplesheet' params(params)

/*
 * Include nf-core modules
 */
include { FASTQC } from './modules/nf-core/fastqc' params(params)
include { MULTIQC } from './modules/nf-core/multiqc' params(params)

/*
 * Run the workflow
 */
workflow {

    CHECK_SAMPLESHEET(ch_input)
        .splitCsv(header:true, sep:',')
        .map { check_samplesheet_paths(it) }
        .set { ch_raw_reads }

    FASTQC(ch_raw_reads)

    OUTPUT_DOCUMENTATION(
        ch_output_docs,
        ch_output_docs_images)

    GET_SOFTWARE_VERSIONS()

    MULTIQC(
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        FASTQC.out.collect(),
        GET_SOFTWARE_VERSIONS.out.yml.collect(),
        ch_workflow_summary)
}

/*
 * Send completion email
 */
workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary, run_name, baseDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}
