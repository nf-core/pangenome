/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input FASTA not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PGGB        } from '../subworkflows/local/pggb'
include { COMMUNITY   } from '../subworkflows/local/community'
include { ODGI_QC     } from '../subworkflows/local/odgi_qc'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { ODGI_SQUEEZE                } from '../modules/nf-core/odgi/squeeze/main'
include { ODGI_VIEW                   } from '../modules/nf-core/odgi/view/main'

//
// MODULE: Locally generated modules
//
include { VG_DECONSTRUCT              } from '../modules/local/vg_deconstruct/main'
include { MULTIQC_COMMUNITY           } from '../modules/local/multiqc_community/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PANGENOME {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    if (params.communities) {
        COMMUNITY (
            INPUT_CHECK.out.fasta,
            INPUT_CHECK.out.fai,
            INPUT_CHECK.out.gzi
        )
        ch_versions = ch_versions.mix(COMMUNITY.out.versions)
        ch_community_join = COMMUNITY.out.fasta_gz.join(COMMUNITY.out.gzi).join(COMMUNITY.out.fai)
        PGGB(
            ch_community_join.map{meta, fasta, gzi, fai -> [ meta, fasta ]},
            ch_community_join.map{meta, fasta, gzi, fai -> [ meta, fai ]},
            ch_community_join.map{meta, fasta, gzi, fai -> [ meta, gzi ]},
            true
        )
        ch_versions = ch_versions.mix(PGGB.out.versions)
        ch_squeeze_in = PGGB.out.og.map{meta, og -> [ [ id:meta.id.replaceFirst(".community.*", "") ], og ]}.groupTuple(by:0)
        ODGI_SQUEEZE(ch_squeeze_in)
        ch_odgi_qc_in = ODGI_SQUEEZE.out.graph.map{meta, ogs -> [ [ id:meta.id + ".squeeze" ], ogs ]}
        ODGI_QC(Channel.empty(), ch_odgi_qc_in, false)
        ch_versions = ch_versions.mix(ODGI_QC.out.versions)
        ch_versions = ch_versions.mix(ODGI_SQUEEZE.out.versions)
        ODGI_VIEW(ch_odgi_qc_in)
        ch_versions = ch_versions.mix(ODGI_VIEW.out.versions)
        ch_multiqc_in = PGGB.out.qc.map{meta, seqwish, gfaffix, viz, viz_pos, viz_depth, viz_inv, viz_O, viz_uncalled, draw -> [ meta, [ seqwish, gfaffix, viz, viz_pos, viz_depth, viz_inv, viz_O, viz_uncalled, draw ] ]}
        MULTIQC_COMMUNITY(ch_multiqc_in,
                            ch_multiqc_config.toList(),
                            ch_multiqc_custom_config.toList(),
                            ch_multiqc_logo.toList())
    } else {
        PGGB (
            INPUT_CHECK.out.fasta,
            INPUT_CHECK.out.fai,
            INPUT_CHECK.out.gzi,
            false
        )
        ch_versions = ch_versions.mix(PGGB.out.versions)
    }

    if (params.vcf_spec != null) {
        ch_vcf_spec = Channel.from(params.vcf_spec).splitCsv().flatten()
        if (!params.communities) {
            VG_DECONSTRUCT(PGGB.out.gfa.combine(ch_vcf_spec))
            ch_versions = ch_versions.mix(VG_DECONSTRUCT.out.versions)
        } else {
            VG_DECONSTRUCT(ODGI_VIEW.out.gfa.combine(ch_vcf_spec))
            ch_versions = ch_versions.mix(VG_DECONSTRUCT.out.versions)
        }

    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPangenome.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPangenome.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    if (!params.communities) {
        if (!params.wfmash_only) {
            ch_multiqc_files = ch_multiqc_files.mix(PGGB.out.qc.map{return it[1..8]})
        }
    } else {
        ch_multiqc_files = ch_multiqc_files.mix(ODGI_QC.out.qc.map{return it[1..8]})
    }

    if (params.vcf_spec != null) {
        ch_multiqc_files = ch_multiqc_files.mix(VG_DECONSTRUCT.out.stats)
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
