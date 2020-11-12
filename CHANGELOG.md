# nf-core/pangenome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/pangenome, created with the [nf-core](https://nf-co.re/) template.

### `Changed`

- `nextflow.config`: Set default container to `ghcr.io/pangenome/pggb:latest`
- `README.md`: Describe how developers can test the current state of the pipeline.
- `test.config`: Point to valid test data and update the default container to `ghcr.io/pangenome/pggb:latest`.

### `Added`

- `main.nf`: Copied the DSL2 pggb workflow from https://github.com/pangenome/pggb/pull/20 into this one. Currently we are ignoring the `nf-core` way to handle input parameters, software versions, etc. This will change in the future.

### `Fixed`

### `Dependencies`

### `Deprecated`
