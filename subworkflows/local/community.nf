//
// Run the Community pipeline
//

include { WFMASH as WFMASH_MAP_COMMUNITY  } from '../../modules/nf-core/wfmash/main'

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

    WFMASH_MAP_COMMUNITY(fasta,
                        query_self,
                        gzi,
                        fai,
                        [],
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

    EXTRACT_COMMUNITIES.out.community_fasta.view()

/*
    WFMASH_MAP(ch_fasta, fai_path, gzi_path, wfmash_prefix)
    paf2Net(WFMASH_MAP.out)
    net2Communities(paf2Net.out)
    extractCommunities(fasta.combine(net2Communities.out.flatten()))
    ch_bgzip_extract_communities = bgzip(fasta.combine(extractCommunities.out.flatten()))
*/

    emit:
    communities = ch_communities
    versions = ch_versions   // channel: [ versions.yml ]
}
