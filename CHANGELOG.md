# nf-core/pangenome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.3 - marsupial

This release brings several template updates up to 3.2.0 and a number of tool updates.

- Templete updates to 3.2.0 @heuermh @subwaystation
- Bump odgi to 0.9.0 @heuermh @subwaystation
- Bump MultiQC to 1.27 @heuermh @subwaystation
- Bump local MultiQC to 1.27 @subwaystation @heuermh
- Bump vg to 1.62 @subwaystation @heuermh
- Bump smoothxg to 0.8.0 @heuermh @subwaystation
- Refactor authors list of contributors @heuermh @subwaystation @Zethson @AndreaGuarracino @mirpedrol @FriederikeHanssen

## 1.1.2 - canguro

This release reverts the `wfmash` tool version to 0.10.4, because the current releases are unstable, not documented, and I don't understand how to set the parameters properly.
It is currently unclear, when a new _stable_ release will become available. @baozg or @AndreaGuarracino are evaluating.

- Per default, we set the number of mappings in `wfmash` to `--n_haplotypes - 1`. Previously, this was set to `1`.
- Deactivated parameter `wfmash_hg_filter_ani_diff`, because the older `wfmash` version does not support it.

## 1.1.1 - LATÜRNICH

This release fixes some important bugs:

- Per default, we set the number of mappings in `wfmash` to `1`. Previously, this was set to the given number of haplotypes.
- To complement the issue above, there is a new parameter `wfmash_n_mappings` with default `1`.
- `bcftools` in the `VG_DECONSTRUCT` module was updated to the most recent version `1.19` to prevent errors like `corrupted size vs. prev_size`.
- Fixed some problems with the delimiter in `VG_DECONSTRUCT` so that the parameters given by `--vcf-spec` are now applied correctly.
- Updated to nf-core template version 2.13.1

## 1.1.0 - Schmuddlweddr

This release mostly contains updates of the respective modules. Functional updates are the following:

- Added parameter `wfmash_hg_filter_ani_diff` mirroring https://github.com/pangenome/pggb/pull/351: Filter out mappings unlikely to be this Average Nuclotide Identitiy (ANI) less than the best mapping. The default values is `30`.
- `seqwish_min_match_length` new default value is `23` now.
- Fixed a bug with `wfmash_exclude_delim` where `#` was interpreted as a comment instead of a string.
- @heringerp added functionality so that very large integer values can be written in short form: For example `5000` would become `5k`.

## 1.0.0 - Ratzupaltuff

Initial release of nf-core/pangenome, created with the [nf-core](https://nf-co.re/) template.
