//
// Check input FASTA and prepare indices
//

include { TABIX_BGZIP                 } from '../../modules/nf-core/tabix/bgzip/main.nf'
include { SAMTOOLS_FAIDX              } from '../../modules/nf-core/samtools/faidx/main.nf'

workflow INPUT_CHECK {
    take:
    fasta // file: /path/to/sequences.fasta

    main:

    ch_versions = Channel.empty() // we collect all versions here
    ch_fasta = Channel.empty() // final output channel [ val(meta) , [ fasta ] ]

    fai_path = file("${params.input}.fai")
    gzi_path = file("${params.input}.gzi")

    fai = Channel.empty() // we store the .fai index here [ fai ]
    gzi = Channel.empty() // we store the .gzi index here [ gzi ]

    meta_fasta = Channel.empty() // intermediate channel where we build our [ val(meta) , [ fasta ] ]
    fasta_file_name = fasta.getName()

    if (params.input.endsWith(".gz")) {
        meta_fasta = tuple([ id:fasta_file_name ], fasta)
        // TODO We want to check, if the input file was actually compressed with bgzip with the upcoming grabix module.
        // For now we assume it was bgzip. If not WFMASH will complain instantly anyhow.
        if (!fai_path.exists() || !gzi_path.exists()) { // the assumption is that none of these files exist if only one does not exist
            SAMTOOLS_FAIDX(meta_fasta, [[],[]])
            fai = SAMTOOLS_FAIDX.out.fai
            gzi = SAMTOOLS_FAIDX.out.gzi
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        } else {
            fai = Channel.of([ [ id:fasta_file_name ], fai_path ])
            gzi = Channel.of([ [ id:fasta_file_name ], gzi_path ])
        }
        ch_fasta = meta_fasta
    } else {
        if (params.input.endsWith("fa")) {
            fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 3)
        } else {
            if (params.input.endswith("fasta")) {
                fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 6)
            } else { // we assume "fna" here
                fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 4)
            }
        }
        meta_fasta = tuple([ id:fasta_file_name ], fasta)
        TABIX_BGZIP(meta_fasta)
        ch_fasta = TABIX_BGZIP.out.output
        SAMTOOLS_FAIDX(ch_fasta, [[],[]])
        gzi = SAMTOOLS_FAIDX.out.gzi
        fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
    }

    emit:
    fasta = ch_fasta         // channel: [ val(meta), [ fasta ] ]
    fai = fai                // channel: [ val(meta), fasta.fai ]
    gzi = gzi                // channel: [ val(meta), fasta.gzi ]
    versions = ch_versions   // channel: [ versions.yml ]
}
