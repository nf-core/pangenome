//
// Run the default PanGenome Graph Builder (PGGB) pipeline
//

include { WFMASH                    } from '../../modules/nf-core/wfmash/main'
include { SEQWISH_INDUCE as SEQWISH } from '../../modules/nf-core/seqwish/induce/main'
include { SMOOTHXG                  } from '../../modules/nf-core/smoothxg/main'
include { GFAFFIX                   } from '../../modules/nf-core/gfaffix/main'
include { ODGI_BUILD                } from '../../modules/nf-core/odgi/build/main'
include { ODGI_UNCHOP               } from '../../modules/nf-core/odgi/unchop/main'
include { ODGI_SORT                 } from '../../modules/nf-core/odgi/sort/main'
include { ODGI_VIEW                 } from '../../modules/nf-core/odgi/view/main'


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
    SEQWISH(ch_seqwish_input)
    ch_versions = ch_versions.mix(SEQWISH.out.versions)

    SMOOTHXG(SEQWISH.out.gfa)
    ch_versions = ch_versions.mix(SMOOTHXG.out.versions)

    GFAFFIX(SMOOTHXG.out.gfa)
    ch_versions = ch_versions.mix(GFAFFIX.out.versions)

    ch_gfaffix = GFAFFIX.out.gfa.map{meta, gfa -> [ [ id: gfa.baseName ], gfa ]}
    ch_seqwish = SEQWISH.out.gfa.map{meta, gfa -> [ [ id: gfa.baseName ], gfa ]}
    ODGI_BUILD(ch_gfaffix.mix(ch_seqwish))
    ch_versions = ch_versions.mix(ODGI_BUILD.out.versions)

    /// this continues the original gfaffix bash script of pggb
    ch_odgi_build = ODGI_BUILD.out.og.map{meta, gfa ->
        if(gfa.baseName.contains("gfaffix")) {
            return [ meta, gfa ]
        } else {
            return null
        }}.filter{ it != null }
    ODGI_UNCHOP(ch_odgi_build)
    ch_versions = ch_versions.mix(ODGI_UNCHOP.out.versions)
    ODGI_SORT(ODGI_UNCHOP.out.unchopped_graph)
    ch_versions = ch_versions.mix(ODGI_SORT.out.versions)
    ODGI_VIEW(ODGI_SORT.out.sorted_graph)
    ch_versions = ch_versions.mix(ODGI_VIEW.out.versions)

    /// TODO GRAPH_QC subworkflow START

    // ODGI_STATS seqwish + sorted gfaffix

    // ODGI_VIZ modes from sorted gfaffix

    // ODGI_Draw modes from sorted gfaffix
    
    /// TODO GRAPH_QC subworkflow END

    emit:
    // TODO qc channel: [ [ seqwish.og.stats.yaml ], [ gfaffix.og.stats.yaml ] ] -> odgi draw and all odgi viz outputs should be going here
    versions = ch_versions   // channel: [ versions.yml ]
}
