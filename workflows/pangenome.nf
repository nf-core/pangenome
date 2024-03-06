/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                 } from 'plugin/nf-validation'
include { paramsSummaryMultiqc             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'

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

workflow PANGENOME {

    main:

    ch_versions = Channel.empty()
    ch_multiqc_report = Channel.empty()

    //
    // SUBWORKFLOW: Read in FASTA, validate and generate appropriate indices
    //
    INPUT_CHECK (
        file(params.input, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()

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
            true,
            false
        )
        ch_versions = ch_versions.mix(PGGB.out.versions)
        ch_squeeze_in = PGGB.out.og.map{meta, og -> [ [ id:meta.id.replaceFirst(".community.*", "") ], og ]}.groupTuple(by:0)
        ODGI_SQUEEZE(ch_squeeze_in)
        ch_odgi_qc_in = ODGI_SQUEEZE.out.graph.map{meta, ogs -> [ [ id:meta.id + ".squeeze" ], ogs ]}
        ODGI_QC(Channel.empty(), ch_odgi_qc_in, false, true)
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
            false,
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

    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pangenome_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    if (!params.communities) {
        if (!params.wfmash_only) {
            ch_multiqc_files = ch_multiqc_files.mix(PGGB.out.qc.map{return it[1..9]})
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

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
