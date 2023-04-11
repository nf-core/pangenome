//
// Run the default PanGenome Graph Builder (PGGB) pipeline
//

include { WFMASH as WFMASH_MAP_ALIGN      } from '../../modules/nf-core/wfmash/main'
include { WFMASH as WFMASH_MAP            } from '../../modules/nf-core/wfmash/main'
include { WFMASH as WFMASH_ALIGN          } from '../../modules/nf-core/wfmash/main'
include { SEQWISH_INDUCE as SEQWISH       } from '../../modules/nf-core/seqwish/induce/main'
include { SMOOTHXG                        } from '../../modules/nf-core/smoothxg/main'
include { GFAFFIX                         } from '../../modules/nf-core/gfaffix/main'
include { ODGI_BUILD                      } from '../../modules/nf-core/odgi/build/main'
include { ODGI_UNCHOP                     } from '../../modules/nf-core/odgi/unchop/main'
include { ODGI_SORT                       } from '../../modules/nf-core/odgi/sort/main'
include { ODGI_VIEW                       } from '../../modules/nf-core/odgi/view/main'

include { SPLIT_APPROX_MAPPINGS_IN_CHUNKS } from '../../modules/local/split_approx_mappings_in_chunks/main'

workflow PGGB {
    take:
    fasta // file: /path/to/sequences.fasta
    fai   // file: /path/to/sequences.fasta.fai
    gzi   // file: /path/to/sequences.fasta.gzi

    main:

    ch_versions = Channel.empty()  // we collect all versions here
    ch_graph_qc = Channel.empty()  // we collect all graph quality control PNGs and graph statistics here
    ch_odgi_view = Channel.empty() // we collect the final graph in GFA format here

    def query_self = true
    if (params.wfmash_only) {
        if (params.wfmash_chunks == 1) {
            WFMASH_MAP_ALIGN(fasta,
                    query_self,
                    gzi,
                    fai,
                    [],
                    [])
            ch_versions = ch_versions.mix(WFMASH_MAP_ALIGN.out.versions)
        } else {
            WFMASH_MAP(fasta,
                    query_self,
                    gzi,
                    fai,
                    [],
                    [])
            ch_versions = ch_versions.mix(WFMASH_MAP.out.versions)
            SPLIT_APPROX_MAPPINGS_IN_CHUNKS(WFMASH_MAP.out.paf)
            ch_versions = ch_versions.mix(SPLIT_APPROX_MAPPINGS_IN_CHUNKS.out.versions)
            ch_split_approx_mappings_in_chunks = SPLIT_APPROX_MAPPINGS_IN_CHUNKS.out.pafs.map{meta, paf -> [ paf ]}.flatten()
            WFMASH_ALIGN(fasta,
                    query_self,
                    gzi,
                    fai,
                    [],
                    ch_split_approx_mappings_in_chunks)
        }
    } else {
        if (params.seqwish_paf != null) {
            file_fasta = file(params.input)
            ch_combined = Channel.of([ [ id:file_fasta.getName() ], file(params.seqwish_paf), file_fasta ])
            SEQWISH(ch_combined)
        } else {
            if (params.wfmash_chunks == 1) {
                WFMASH_MAP_ALIGN(fasta,
                        query_self,
                        gzi,
                        fai,
                        [],
                        [])
                ch_versions = ch_versions.mix(WFMASH_MAP_ALIGN.out.versions)
                ch_seqwish_input = WFMASH_MAP_ALIGN.out.paf.join(fasta)
                SEQWISH(ch_seqwish_input)
                ch_versions = ch_versions.mix(SEQWISH.out.versions)
            } else {
                WFMASH_MAP(fasta,
                        query_self,
                        gzi,
                        fai,
                        [],
                        [])
                ch_versions = ch_versions.mix(WFMASH_MAP.out.versions)
                SPLIT_APPROX_MAPPINGS_IN_CHUNKS(WFMASH_MAP.out.paf)
                ch_versions = ch_versions.mix(SPLIT_APPROX_MAPPINGS_IN_CHUNKS.out.versions)
                ch_split_approx_mappings_in_chunks = SPLIT_APPROX_MAPPINGS_IN_CHUNKS.out.pafs.map{meta, paf -> [ paf ]}.flatten()
                WFMASH_ALIGN(fasta,
                        query_self,
                        gzi,
                        fai,
                        [],
                        ch_split_approx_mappings_in_chunks)
                SEQWISH(WFMASH_ALIGN.out.paf.groupTuple(by: 0, size: params.wfmash_chunks).join(fasta))
            }
        }

        if (params.skip_smoothxg) {
            GFAFFIX(SEQWISH.out.gfa)
            ch_versions = ch_versions.mix(GFAFFIX.out.versions)
        } else {
            SMOOTHXG(SEQWISH.out.gfa)
            ch_versions = ch_versions.mix(SMOOTHXG.out.versions)
            GFAFFIX(SMOOTHXG.out.gfa)
            ch_versions = ch_versions.mix(GFAFFIX.out.versions)
        }

        ch_gfaffix = GFAFFIX.out.gfa.map{meta, gfa -> [ [ id: gfa.baseName ], gfa ]}
        ch_seqwish = SEQWISH.out.gfa.map{meta, gfa -> [ [ id: gfa.baseName ], gfa ]}
        ODGI_BUILD(ch_gfaffix.mix(ch_seqwish))
        ch_versions = ch_versions.mix(ODGI_BUILD.out.versions)

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
        ch_odgi_view = ODGI_VIEW.out.gfa

        ch_odgi_build_seqwish = ODGI_BUILD.out.og.map{meta, gfa ->
            if(gfa.baseName.contains("seqwish")) {
                return [ meta, gfa ]
            } else {
                return null
            }}.filter{ it != null }
    }

    emit:
    gfa = ch_odgi_view
    seqwish = ch_odgi_build_seqwish
    sorted_graph = ODGI_SORT.out.sorted_graph
    versions = ch_versions   // channel: [ versions.yml ]
}