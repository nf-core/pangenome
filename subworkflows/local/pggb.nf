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
include { ODGI_QC                         } from '../../subworkflows/local/odgi_qc'

workflow PGGB {
    take:
    fasta // file: /path/to/sequences.fasta
    fai   // file: /path/to/sequences.fasta.fai
    gzi   // file: /path/to/sequences.fasta.gzi
    community_mode // val: determines how we will build our meta identifiers of ODGI_QC
    no_seqwish_input // val : determines how we will sort the input for MULTIQC

    main:

    ch_versions = Channel.empty()  // we collect all versions here
    ch_graph_qc = Channel.empty()  // we collect all graph quality control PNGs and graph statistics here
    ch_odgi_view = Channel.empty() // we collect the final graph in GFA format here
    ch_odgi_build_seqwish = Channel.empty()
    ch_sorted_graph = Channel.empty()

    def query_self = true
    if (params.wfmash_only) {
        if (params.wfmash_chunks == 1) {
            ch_wfmash_map_align = fasta.map{meta, fasta -> [ meta, fasta, []]}
            ch_wfmash_map_align = ch_wfmash_map_align.join(gzi).join(fai)
            WFMASH_MAP_ALIGN(ch_wfmash_map_align,
                    query_self,
                    [])
            ch_versions = ch_versions.mix(WFMASH_MAP_ALIGN.out.versions)
        } else {
            ch_wfmash_map = fasta.map{meta, fasta -> [ meta, fasta, [] ]}
            ch_wfmash_map = ch_wfmash_map.join(gzi).join(fai)
            WFMASH_MAP(ch_wfmash_map,
                    query_self,
                    [])
            ch_versions = ch_versions.mix(WFMASH_MAP.out.versions)
            SPLIT_APPROX_MAPPINGS_IN_CHUNKS(WFMASH_MAP.out.paf)
            ch_versions = ch_versions.mix(SPLIT_APPROX_MAPPINGS_IN_CHUNKS.out.versions)
            ch_wfmash_align = fasta.combine(SPLIT_APPROX_MAPPINGS_IN_CHUNKS.out.pafs.transpose(), by:0).combine(gzi, by:0).combine(fai, by:0)
            WFMASH_ALIGN(ch_wfmash_align,
                    query_self,
                    [])
            ch_versions = ch_versions.mix(WFMASH_ALIGN.out.versions)
        }
    } else {
        if (params.seqwish_paf != null) {
            file_fasta = file(params.input)
            ch_combined = Channel.of([ [ id:file_fasta.getName() ], file(params.seqwish_paf), file_fasta ])
            SEQWISH(ch_combined)
            ch_versions = ch_versions.mix(SEQWISH.out.versions)
        } else {
            if (params.wfmash_chunks == 1) {
                ch_wfmash_map_align = fasta.map{meta, fasta -> [ meta, fasta, [] ]}
                ch_wfmash_map_align = ch_wfmash_map_align.join(gzi).join(fai)
                WFMASH_MAP_ALIGN(ch_wfmash_map_align,
                        query_self,
                        [])
                ch_versions = ch_versions.mix(WFMASH_MAP_ALIGN.out.versions)
                ch_seqwish_input = WFMASH_MAP_ALIGN.out.paf.join(fasta)
                SEQWISH(ch_seqwish_input)
                ch_versions = ch_versions.mix(SEQWISH.out.versions)
            } else {
                ch_wfmash_map = fasta.map{meta, fasta -> [ meta, fasta, [] ]}
                ch_wfmash_map = ch_wfmash_map.join(gzi).join(fai)
                WFMASH_MAP(ch_wfmash_map,
                        query_self,
                        [])
                ch_versions = ch_versions.mix(WFMASH_MAP.out.versions)
                SPLIT_APPROX_MAPPINGS_IN_CHUNKS(WFMASH_MAP.out.paf)
                ch_versions = ch_versions.mix(SPLIT_APPROX_MAPPINGS_IN_CHUNKS.out.versions)
                ch_wfmash_align = fasta.combine(SPLIT_APPROX_MAPPINGS_IN_CHUNKS.out.pafs.transpose(), by:0).combine(gzi, by:0).combine(fai, by:0)
                WFMASH_ALIGN(ch_wfmash_align,
                        query_self,
                        [])
                ch_versions = ch_versions.mix(WFMASH_ALIGN.out.versions)
                SEQWISH(WFMASH_ALIGN.out.paf.groupTuple(by: 0, size: params.wfmash_chunks).join(fasta))
                ch_versions = ch_versions.mix(SEQWISH.out.versions)
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

        ch_sorted_graph = ODGI_SORT.out.sorted_graph

        ODGI_QC(ch_odgi_build_seqwish, ODGI_SORT.out.sorted_graph, community_mode, no_seqwish_input)
        ch_graph_qc = ODGI_QC.out.qc
        ch_versions = ch_versions.mix(ODGI_QC.out.versions)
    }

    emit:
    gfa = ch_odgi_view
    qc = ch_graph_qc
    og = ch_sorted_graph
    versions = ch_versions   // channel: [ versions.yml ]
}
