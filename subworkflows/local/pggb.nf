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
include { ODGI_STATS                      } from '../../modules/nf-core/odgi/stats/main'
include { ODGI_VIZ as ODGI_VIZ_COLOR      } from '../../modules/nf-core/odgi/viz/main'
include { ODGI_VIZ as ODGI_VIZ_POS        } from '../../modules/nf-core/odgi/viz/main'
include { ODGI_VIZ as ODGI_VIZ_DEPTH      } from '../../modules/nf-core/odgi/viz/main'
include { ODGI_VIZ as ODGI_VIZ_INV        } from '../../modules/nf-core/odgi/viz/main'
include { ODGI_VIZ as ODGI_VIZ_UNCALLED   } from '../../modules/nf-core/odgi/viz/main'
include { ODGI_VIZ as ODGI_VIZ_COMPR      } from '../../modules/nf-core/odgi/viz/main'
include { ODGI_LAYOUT                     } from '../../modules/nf-core/odgi/layout/main'
include { ODGI_DRAW as ODGI_DRAW_MULTIQC  } from '../../modules/nf-core/odgi/draw/main'
include { ODGI_DRAW as ODGI_DRAW_HEIGHT   } from '../../modules/nf-core/odgi/draw/main'

include { SPLIT_APPROX_MAPPINGS_IN_CHUNKS } from '../../modules/local/split_approx_mappings_in_chunks/main'

workflow PGGB {
    take:
    fasta // file: /path/to/sequences.fasta
    fai   // file: /path/to/sequences.fasta.fai
    gzi   // file: /path/to/sequences.fasta.gzi

    main:

    ch_versions = Channel.empty() // we collect all versions here
    ch_graph_qc = Channel.empty() // we collect all graph quality control PNGs here

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

        ch_odgi_build_seqwish = ODGI_BUILD.out.og.map{meta, gfa ->
            if(gfa.baseName.contains("seqwish")) {
                return [ meta, gfa ]
            } else {
                return null
            }}.filter{ it != null }
        ODGI_STATS(ch_odgi_build_seqwish.mix(ODGI_SORT.out.sorted_graph))
        ch_versions = ch_versions.mix(ODGI_STATS.out.versions)
        // prepare all inputs
        ch_odgi_viz = ODGI_SORT.out.sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz" ], gfa ]}
        ODGI_VIZ_COLOR(ch_odgi_viz)
        ch_versions = ch_versions.mix(ODGI_VIZ_COLOR.out.versions)
        ch_odgi_viz_pos = ODGI_SORT.out.sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_pos" ], gfa ]}
        ODGI_VIZ_POS(ch_odgi_viz_pos)
        ch_versions = ch_versions.mix(ODGI_VIZ_POS.out.versions)
        ch_odgi_viz_depth = ODGI_SORT.out.sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_depth" ], gfa ]}
        ODGI_VIZ_DEPTH(ch_odgi_viz_depth)
        ch_versions = ch_versions.mix(ODGI_VIZ_DEPTH.out.versions)
        ch_odgi_viz_inv = ODGI_SORT.out.sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_inv" ], gfa ]}
        ODGI_VIZ_INV(ch_odgi_viz_inv)
        ch_versions = ch_versions.mix(ODGI_VIZ_INV.out.versions)
        ch_odgi_viz_compr = ODGI_SORT.out.sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_O" ], gfa ]}
        ODGI_VIZ_COMPR(ch_odgi_viz_compr)
        ch_versions = ch_versions.mix(ODGI_VIZ_COMPR.out.versions)
        ch_odgi_viz_uncalled = ODGI_SORT.out.sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_uncalled" ], gfa ]}
        ODGI_VIZ_UNCALLED(ch_odgi_viz_uncalled)
        ch_versions = ch_versions.mix(ODGI_VIZ_UNCALLED.out.versions)

        ODGI_LAYOUT(ODGI_SORT.out.sorted_graph)
        ch_versions = ch_versions.mix(ODGI_LAYOUT.out.versions)
        ch_odgi_layout = ODGI_LAYOUT.out.lay.map{meta, lay -> [ lay ]}
        ODGI_DRAW_MULTIQC(ODGI_SORT.out.sorted_graph, ch_odgi_layout)
        ch_versions = ch_versions.mix(ODGI_DRAW_MULTIQC.out.versions)
        ODGI_DRAW_HEIGHT(ODGI_SORT.out.sorted_graph, ch_odgi_layout)
        ch_versions = ch_versions.mix(ODGI_DRAW_HEIGHT.out.versions)

        /// TODO GRAPH_QC subworkflow END

        ch_stats = ODGI_STATS.out.yaml.map{meta, yaml -> return yaml}.collect()

        ch_viz = ODGI_VIZ_COLOR.out.png.map{meta, png -> return png}.collect()
        ch_viz_pos = ODGI_VIZ_POS.out.png.map{meta, png -> return png}.collect()
        ch_viz_depth = ODGI_VIZ_DEPTH.out.png.map{meta, png -> return png}.collect()
        ch_viz_inv = ODGI_VIZ_INV.out.png.map{meta, png -> return png}.collect()
        ch_viz_compr = ODGI_VIZ_COMPR.out.png.map{meta, png -> return png}.collect()
        ch_viz_uncalled = ODGI_VIZ_UNCALLED.out.png.map{meta, png -> return png}.collect()

        ch_draw = ODGI_DRAW_MULTIQC.out.png.map{meta, png -> return png}.collect()

        ch_graph_qc = ch_stats.mix(ch_viz, ch_viz_pos, ch_viz_depth, ch_viz_inv, ch_viz_compr, ch_viz_uncalled, ch_draw)

    }

    emit:
    gfa = ODGI_VIEW.out.gfa
    qc = ch_graph_qc // [ seqwish.og.stats.yaml , gfaffix.og.stats.yaml, odgi_viz_multiqc.png, odgi_viz_pos_multiqc.png, odgi_viz_depth_multiqc.png, odgi_viz_inv_multiqc.png, odgi_viz_compr_multiqc.png, odgi_viz_uncalled_multiqc.png, odgi_draw_multiqc.png ]
    versions = ch_versions   // channel: [ versions.yml ]
}
