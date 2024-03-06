//
// Subworkflow with functionality specific to the nf-core/pangenome pipeline
//

import groovy.json.JsonSlurper

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input <BGZIPPED_FASTA> --n_haplotypes <NUM_HAPS_IN_FASTA> --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same strandedness
    def strandedness_ok = metas.collect{ it.strandedness }.unique().size == 1
    if (!strandedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must have the same strandedness!: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    // Check mandatory parameters
    if (!params.input) {
        exit 1, 'Input FASTA not specified!'
    }
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "NF-CORE (Ewels et al. 2020)",
            "Nextflow (Di Tommaso et al. 2017)",
            "GFAFFIX (Liao et al. 2023)",
            "MultiQC (Ewels et al. 2016)",
            "NET2COMMUNITIES (Traag et al. 2019)",
            "ODGI (Guarracino, Heumos et al. 2022; Heumos, Guarracino et al. 2023)",
            "PGGB (Garrison, Guarracino et al. 2023; Guarracino et al. 2023; Liao et al. 2023)",
            "SAMTOOLS (Li et al. 2009)",
            "SEQWISH (Garrison and Guarracino 2023)",
            "SMOOTHXG (Garrison, Guarracino et al. 2023)",
            "VCFLIB (Garrison et al. 2022)",
            "VG (Garrison et al. 2018)",
            "WFMASH (Guarracino et al. 2024)",
            "Anaconda (Anaconda Software Distribution 2016)",
            "Bioconda (Grüning et al. 2018)",
            "BioContainers (da Veiga Leprevost et al. 2017)",
            "Docker (Merkel 2014)",
            "Singularity (Kurtzer et al. 2017)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: https://10.1038/s41587-020-0439-x. PubMed PMID: 32055031.</li>",
            "<li>Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: https://10.1038/nbt.3820. PubMed PMID: 28398311.</li>",
            "<li>Liao W, Asri M, Ebler J, Doerr D, Haukness M, Hickey G, Lu S, Lucas J K, Monlong J, Abel H J, Buonaiuto S, Chang X H, Cheng H, Chu J, Colonna V, Eizenga J M, Feng X, Fischer C, Fulton R S, Garg S, Groza C, Guarracino A, Harvey W T, Heumos S, Howe K, Jain M, Lu T, Markello C, Martin F J, Mitchell M W, Munson K M, Mwaniki M N, Novak A M, Olsen H E, Pesout T, Porubsky D, Prins P, Sibbesen J A, Tomlinson C, Villani F, Vollger M R, Human Pangenome Reference Consortium, Bourque G, Chaisson M J P, Flicek P, Phillippy A M, Zook J M, Eichler E E,Haussler D, Jarvis E D, Miga K H, Wang T, Garrison E, Marschall T, Hall I, Li H, Paten B. A Draft Human Pangenome Reference. Nature 617, 312–324 (2023). https://doi.org/10.1038/s41586-023-05896-x.</li>",
            "<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: https://10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.</li>",
            "<li>Traag, VA, Waltman, L & van Eck, NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep 9, 5233 (2019). https://doi.org/10.1038/s41598-019-41695-z.</li>",
            "<li>Guarracino A, Heumos S, Nahnsen S, Prins P, Garrison E. ODGI: understanding pangenome graphs. Bioinformatics. Volume 38. Issue 13. July 2022. Pages 3319–3326. https://doi.org/10.1093/bioinformatics/btac308.</li>",
            "<li>Heumos S, Guarracino A, Schmelzle J N M, Li J, Zhang Z, Nahnsen S, Prins P, Garrison E. Pangenome graph layout by Path-Guided Stochastic Gradient Descent. bioRxiv. https://www.biorxiv.org/content/10.1101/2023.09.22.558964v1.</li>",
            "<li>Garrison E, Guarracino A, Heumos S, Villani F, Bao Z, Tattini L, Hagmann J, Vorbrugg S, Marco-Sola S, Kubica S, Ashbrook D G, Thorell K, Rusholme-Pilcher R L, Liti G, Rudbeck E, Nahnsen S, Yang Z, Moses M N, Nobrega F L, Wu Y, Chen H, de Ligt J, Sudmant P H, Soranzo N, Colonna V, Williams R W, Prins P. Building pangenome graphs. bioRxiv. 2023.04.05.535718. doi: https://doi.org/10.1101/2023.04.05.535718.</li>",
            "<li>Guarracino A, Buonaiuto S, de Lima L G, Potapova T, Rhie A, Koren S, Rubinstein B, Fischer C, Human Pangenome Reference Consortium, Gerton J L, Phillippy A M, Colonna V, Garrison E. Recombination between heterologous human acrocentric chromosomes. Nature 617, 335–343 (2023). https://doi.org/10.1038/s41586-023-05976-y.</li>",
            "<li>Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: https://10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PMID: 19505943; PMCID: PMC2723002.</li>",
            "<li>Garrison E, Guarracino A. Unbiased pangenome graphs. Bioinformatics. 2023 Jan 1;39(1):btac743. doi: https://10.1093/bioinformatics/btac743. PMID: 36448683; PMCID: PMC9805579.</li>",
            "<li>Garrison E, Kronenberg ZN, Dawson ET, Pedersen BS, Prins P (2022) A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar. PLOS Computational Biology 18(5): e1009123. https://doi.org/10.1371/journal.pcbi.1009123.</li>",
            "<li>Garrison E et al. Variation graph toolkit improves read mapping by representing genetic variation in the reference. Nature biotechnology vol. 36,9 (2018): 875-879. doi: https://10.1038/nbt.4227.</li>",
            "<li>Guarracino A, Mwaniki N, Marco-Sola S, Garrison E. Wfmash: whole-chromosome pairwise alignment using the hierarchical wavefront algorithm. 2024. https://github.com/waveygang/wfmash.</li>",
            "<li>Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.</li>",
            "<li>Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: https://10.1038/s41592-018-0046-7. PubMed PMID: 29967506.</li>",
            "<li>da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: https://10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.</li>",
            "<li>Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: https://10.5555/2600239.2600241.</li>",
            "<li>Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: https://10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()
    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
