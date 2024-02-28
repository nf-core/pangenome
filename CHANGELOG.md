# nf-core/pangenome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - Schmuddlweddr

This release mostly contains updates of the respective modules. Functional updates are the following:

- Added parameter `wfmash_hg_filter_ani_diff` mirroring https://github.com/pangenome/pggb/pull/351: Filter out mappings unlikely to be this Average Nuclotide Identitiy (ANI) less than the best mapping. The default values is `30`.
- `seqwish_min_match_length` new default value is `23` now.
- Fixed a bug with `wfmash_exclude_delim` where `#` was interpreted as a comment instead of a string.
- @heringerp added functionality so that very large integer values can be written in short form: For example `5000` would become `5k`.

## 1.0.0 - Ratzupaltuff

Initial release of nf-core/pangenome, created with the [nf-core](https://nf-co.re/) template.
