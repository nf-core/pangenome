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
        meta = [ id:fasta_file_name ]
        meta_fasta = tuple(meta, fasta)
        // TODO We want to check, if the input file was actually compressed with bgzip with the upcoming grabix module.
        // For now we assume it was bgzip. If not wfmash will complain instantly anyhow.
        if (!fai_path.exists() || !gzi_path.exists()) { // the assumption is that none of the files exist if only one does not exist
            SAMTOOLS_FAIDX(meta_fasta)
            fai = SAMTOOLS_FAIDX.out.fai.collect{it[1]}
            gzi = SAMTOOLS_FAIDX.out.gzi.collect{it[1]}
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        } else {
            fai = Channel.fromPath("${params.input}.fai").collect()
            gzi = Channel.fromPath("${params.input}.gzi").collect()
        }
        ch_fasta = meta_fasta
    } else {
        if (params.input.endsWith("fa")) {
            // SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].baseName], it] })
            // https://stackoverflow.com/questions/20954779/regular-expression-to-get-last-3-characters-of-a-string
            // ... id:get_last_three_chars(it[0]) ...
            fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 3)
        } else {
            if (params.input.endswith("fasta")) {
                fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 6)
            } else { // we assume "fna" here
                fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 4)
            }
        }
        meta = [ id:fasta_file_name ]
        meta_fasta = tuple(meta, fasta)
        TABIX_BGZIP(meta_fasta)
        ch_fasta = TABIX_BGZIP.out.output
        SAMTOOLS_FAIDX(ch_fasta)
        gzi = SAMTOOLS_FAIDX.out.gzi.collect{it[1]}
        fai = SAMTOOLS_FAIDX.out.fai.collect{it[1]}
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
    }

    emit:
    fasta = ch_fasta         // channel: [ val(meta), [ fasta ] ]
    fai = fai                // channel: [ fasta.fai ]
    gzi = gzi                // channel: [ fasta.gzi ]
    versions = ch_versions   // channel: [ versions.yml ]
}
