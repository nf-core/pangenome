//
// Run the default PanGenome Graph Builder (PGGB) pipeline
//

include { WFMASH                    } from '../../modules/nf-core/wfmash/main'
include { SEQWISH_INDUCE as SEQWISH } from '../../modules/nf-core/seqwish/induce/main'
include { SMOOTHXG                  } from '../../modules/nf-core/smoothxg/main'
include { GFAFFIX                   } from '../../modules/nf-core/gfaffix/main'

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

    ch_seqwish_input = WFMASH.out.paf.join(fasta) // TODO I want to have meta.id = meta.id + ".seqwish"
    SEQWISH(ch_seqwish_input) // tuple val(meta), path("*.gfa"), emit: gfa
    ch_versions = ch_versions.mix(SEQWISH.out.versions)

    SMOOTHXG(SEQWISH.out.gfa)
    ch_versions = ch_versions.mix(SMOOTHXG.out.versions)

    GFAFFIX(SMOOTHXG.out.gfa)
    ch_versions = ch_versions.mix(GFAFFIX.out.versions)

    
    // TODO ODGI_BUILD of SEQWISH, GFAFFIX

    // TODO ODGI_UNCHOP

    // TODO ODGI_SORT

    // TODO ODGI_VIEW

    emit:
    // TODO gfaffix channel: [ [meta.id], [ graph.og ] ] -> we actually take the ODGI_SORT output!
    // TODO seqwish channel: [ [meta.id], [ graph.og ] ]
    versions = ch_versions   // channel: [ versions.yml ]
}
