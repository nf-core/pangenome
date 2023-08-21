//
// Run the Community pipeline
//

include { WFMASH as WFMASH_MAP_COMMUNITY  } from '../../modules/nf-core/wfmash/main'
include { TABIX_BGZIP                     } from '../../modules/nf-core/tabix/bgzip/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main.nf'

include { PAF2NET             } from '../../modules/local/paf2net/main'
include { NET2COMMUNITIES     } from '../../modules/local/net2communities/main'
include { EXTRACT_COMMUNITIES } from '../../modules/local/extract_communities/main'

workflow COMMUNITY {
    take:
    fasta // file: /path/to/sequences.fasta
    fai   // file: /path/to/sequences.fasta.fai
    gzi   // file: /path/to/sequences.fasta.gzi

    main:

    ch_versions = Channel.empty() // we collect all versions here
    ch_communities = Channel.empty()

    def query_self = true

    ch_wfmash_map = fasta.map{meta, fasta -> [ meta, fasta, [] ]}
    ch_wfmash_map = ch_wfmash_map.join(gzi).join(fai)
    WFMASH_MAP_COMMUNITY(ch_wfmash_map,
                        query_self,
                        [])
    ch_versions = ch_versions.mix(WFMASH_MAP_COMMUNITY.out.versions)

    PAF2NET(WFMASH_MAP_COMMUNITY.out.paf)
    ch_versions = ch_versions.mix(PAF2NET.out.versions)

    NET2COMMUNITIES(PAF2NET.out.txts)
    ch_versions = ch_versions.mix(NET2COMMUNITIES.out.versions)

    ch_txt_communities = fasta.combine(NET2COMMUNITIES.out.communities.flatten())
    ch_txt_communities = ch_txt_communities.map{meta, fasta, community -> [ [ id: community.baseName.split("//.")[-1] ], fasta, community ]}

    EXTRACT_COMMUNITIES(ch_txt_communities)
    ch_versions = ch_versions.mix(EXTRACT_COMMUNITIES.out.versions)

    TABIX_BGZIP(EXTRACT_COMMUNITIES.out.community_fasta)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)

    SAMTOOLS_FAIDX(TABIX_BGZIP.out.output, [[],[]])
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    fasta_gz = TABIX_BGZIP.out.output         // channel: [ val(meta), [ fasta.gz ] ]
    gzi = SAMTOOLS_FAIDX.out.gzi
    fai = SAMTOOLS_FAIDX.out.fai
    versions = ch_versions   // channel: [ versions.yml ]
}
