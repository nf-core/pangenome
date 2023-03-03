//
// Run the default PanGenome Graph Builder (PGGB) pipeline
//

include { WFMASH                    } from '../../modules/nf-core/wfmash/main'
include { SEQWISH_INDUCE as SEQWISH } from '../../modules/nf-core/seqwish/induce/main'
include { SMOOTHXG                  } from '../../modules/nf-core/smoothxg/main'

workflow PGGB {
    take:
    fasta // file: /path/to/sequences.fasta
    fai   // file: /path/to/sequences.fasta.fai
    gzi   // file: /path/to/sequences.fasta.gzi

    main:

    ch_versions = Channel.empty() // we collect all versions here

    def query_self = true
    WFMASH(fasta,
            query_self,
            gzi,
            fai,
            [],
            [])
    ch_versions = ch_versions.mix(WFMASH.out.versions)

    ch_seqwish_input = WFMASH.out.paf.join(fasta)
    SEQWISH(ch_seqwish_input) // tuple val(meta), path("*.gfa"), emit: gfa
    ch_versions = ch_versions.mix(SEQWISH.out.versions)

    SMOOTHXG(SEQWISH.out.gfa)
    ch_versions = ch_versions.mix(SMOOTHXG.out.versions)

    emit:
    versions = ch_versions   // channel: [ versions.yml ]
}
