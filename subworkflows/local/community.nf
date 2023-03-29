//
// Run the Community pipeline
//

include { WFMASH as WFMASH_MAP_COMMUNITY  } from '../../modules/nf-core/wfmash/main'

include { PAF_2_NET } from '../../modules/local/paf2net/main'

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

    PAF_2_NET(WFMASH_MAP_COMMUNITY.out.paf)
    ch_versions = ch_versions.mix(PAF_2_NET.out.versions)

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
