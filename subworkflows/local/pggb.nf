//
// Run the default PanGenome Graph Builder (PGGB) pipeline
//

include { WFMASH                    } from '../../modules/nf-core/wfmash/main'
include { SEQWISH_INDUCE as SEQWISH } from '../../modules/nf-core/seqwish/induce/main'
include { SMOOTHXG                  } from '../../modules/nf-core/smoothxg/main'
include { GFAFFIX                   } from '../../modules/nf-core/gfaffix/main'
include { ODGI_BUILD                } from '../../modules/nf-core/odgi/build/main'

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

    ch_gfaffix = GFAFFIX.out.gfa.map{meta, gfa -> [ [ id: gfa.baseName ], gfa ]}
    ch_seqwish = SEQWISH.out.gfa.map{meta, gfa -> [ [ id: gfa.baseName ], gfa ]}
    ODGI_BUILD(ch_gfaffix.mix(ch_seqwish))
    ch_versions = ch_versions.mix(ODGI_BUILD.out.versions)

    // TODO ODGI_UNCHOP
    ch_odgi_build = ODGI_BUILD.out.og.map{meta, gfa ->
        if(gfa.baseName.contains("gfaffix")) {
            return [ meta, gfa ]
        } else {
            return null
        }}.filter{ it != null }
    // TODO ODGI_SORT

    // TODO ODGI_VIEW

    emit:
    // TODO gfaffix channel: [ [meta.id], [ graph.og ] ] -> we actually take the ODGI_SORT output!
    // TODO seqwish channel: [ [meta.id], [ graph.og ] ]
    versions = ch_versions   // channel: [ versions.yml ]
}
