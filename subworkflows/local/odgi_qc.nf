//
// Run the default ODGI quality control subworkflow
//

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

workflow ODGI_QC {
    take:
    seqwish          // file: /path/to/seqwish.og
    sorted_graph     // file: /path/to/sorted_graph.og
    community_mode   // val : determines how to handly meta identifiers
    no_seqwish_input // val : determines how we will sort the input for MULTIQC

    main:

    ch_versions = Channel.empty() // we collect all versions here
    ch_graph_qc = Channel.empty() // we collect all graph quality control PNGs here

    ODGI_STATS(seqwish.mix(sorted_graph))
    ch_versions = ch_versions.mix(ODGI_STATS.out.versions)
    // prepare all inputs
    ch_odgi_viz = sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz" ], gfa ]}
    ODGI_VIZ_COLOR(ch_odgi_viz)
    ch_versions = ch_versions.mix(ODGI_VIZ_COLOR.out.versions)
    ch_odgi_viz_pos = sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_pos" ], gfa ]}
    ODGI_VIZ_POS(ch_odgi_viz_pos)
    ch_versions = ch_versions.mix(ODGI_VIZ_POS.out.versions)
    ch_odgi_viz_depth = sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_depth" ], gfa ]}
    ODGI_VIZ_DEPTH(ch_odgi_viz_depth)
    ch_versions = ch_versions.mix(ODGI_VIZ_DEPTH.out.versions)
    ch_odgi_viz_inv = sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_inv" ], gfa ]}
    ODGI_VIZ_INV(ch_odgi_viz_inv)
    ch_versions = ch_versions.mix(ODGI_VIZ_INV.out.versions)
    ch_odgi_viz_compr = sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_O" ], gfa ]}
    ODGI_VIZ_COMPR(ch_odgi_viz_compr)
    ch_versions = ch_versions.mix(ODGI_VIZ_COMPR.out.versions)
    ch_odgi_viz_uncalled = sorted_graph.map{meta, gfa -> [ [ id: meta.id + ".viz_uncalled" ], gfa ]}
    ODGI_VIZ_UNCALLED(ch_odgi_viz_uncalled)
    ch_versions = ch_versions.mix(ODGI_VIZ_UNCALLED.out.versions)

    ODGI_LAYOUT(sorted_graph)
    ch_versions = ch_versions.mix(ODGI_LAYOUT.out.versions)
    ODGI_DRAW_MULTIQC(sorted_graph.join(ODGI_LAYOUT.out.lay))
    ch_versions = ch_versions.mix(ODGI_DRAW_MULTIQC.out.versions)
    ODGI_DRAW_HEIGHT(sorted_graph.join(ODGI_LAYOUT.out.lay))
    ch_versions = ch_versions.mix(ODGI_DRAW_HEIGHT.out.versions)

    if (community_mode) {
        ch_graph_qc = ODGI_STATS.out.yaml.map{meta, ymls -> [ [ id: meta.id.replaceFirst(".seqwish", "").replaceFirst(".gfaffix", "") ], ymls ]}.groupTuple(by:0, size:2).map{meta, ymls -> [ meta, ymls[0], ymls[1] ]}
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_COLOR.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".gfaffix.*", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_POS.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".gfaffix.*", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_DEPTH.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".gfaffix.*", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_INV.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".gfaffix.*", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_COMPR.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".gfaffix.*", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_UNCALLED.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".gfaffix.*", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_DRAW_MULTIQC.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".gfaffix.*", "") ], png ]})
    } else {
        //
        if (no_seqwish_input) {
            ch_graph_qc = ODGI_STATS.out.yaml
        } else {
            ch_graph_qc = ODGI_STATS.out.yaml.collect().map{seqwish_id, seqwish_qc, gfaffix_id, gfaffix_qc -> [ gfaffix_id, gfaffix_qc, seqwish_qc ]}
        }
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_COLOR.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".viz", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_POS.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".viz_pos", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_DEPTH.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".viz_depth", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_INV.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".viz_inv", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_COMPR.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".viz_O", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_VIZ_UNCALLED.out.png.map{meta, png -> [ [ id: meta.id.replaceFirst(".viz_uncalled", "") ], png ]})
        ch_graph_qc = ch_graph_qc.join(ODGI_DRAW_MULTIQC.out.png)
    }

    emit:
    qc = ch_graph_qc // [ seqwish.og.stats.yaml , gfaffix.og.stats.yaml, odgi_viz_multiqc.png, odgi_viz_pos_multiqc.png, odgi_viz_depth_multiqc.png, odgi_viz_inv_multiqc.png, odgi_viz_compr_multiqc.png, odgi_viz_uncalled_multiqc.png, odgi_draw_multiqc.png ]
    versions = ch_versions   // channel: [ versions.yml ]
}
