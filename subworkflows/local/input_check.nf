//
// Check input FASTA and prepare indices
//

include { TABIX_BGZIP                } from '../../modules/nf-core/tabix/bgzip/main.nf'
include { SAMTOOLS_FAIDX             } from '../../modules/nf-core/samtools/faidx/main.nf'

workflow INPUT_CHECK {
    take:
    fasta // file: /path/to/sequences.fasta

    main:

    fai_path = file("${params.input}.fai")
    gzi_path = file("${params.input}.gzi")

    def fasta_file_name = fasta.getName()
    // TODO this only applies when we end with *.fa!
    fasta_file_name = fasta_file_name.substring(0, fasta_file_name.length() - 3)
    meta = [ id:fasta_file_name ]
    ch_fasta = tuple(meta, fasta)

    TABIX_BGZIP(ch_fasta)
    SAMTOOLS_FAIDX(TABIX_BGZIP.out.output)

    // TODO wfmash testing
    /*
    input:
    tuple val(meta), path(fasta_gz)
    val(query_self)
    path(gzi)
    path(fai)
    path(fasta_query_list)
    path(paf)
    */

/*
    if (params.input.endsWith(".gz")) {
        if (!fai_path.exists() || !gzi_path.exists()) { // the assumption is that none of the files exist if only one does not exist
        // samtoolsFaidx(fasta)
        // fai = samtoolsFaidx.out.samtools_fai.collect()
        // gzi = samtoolsFaidx.out.samtools_gzi.collect()
        } else {
        fai = channel.fromPath("${params.input}.fai").collect()
        gzi = channel.fromPath("${params.input}.gzi").collect()
        }
    } else {
        // val(meta), path(input)
        TABIX_BGZIP(ch_fasta)

    }
*/
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // TODO combine output so that we have

    emit:
    fasta                  // channel: [ val(meta), [ reads ] ]
    versions = ch_versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

