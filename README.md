# ![nf-core/pangenome](docs/images/nf-core-pangenome_logo.png)

**The pangenome graph construction pipeline renders a collection of sequences into a pangenome graph. Its goal is to build a graph that is locally directed and acyclic while preserving large-scale variation. Maintaining local linearity is important for interpretation, visualization, mapping, comparative genomics, and reuse of pangenome graphs**.

> As a first step, the aim is to create a Nextflow DSL2 version of the pangenome graph builder [`pggb`](https://github.com/pangenome/pggb) pipeline.

> This pipeline is still under heavy development! The code is likely to change and the documentation is imperfectly realized. An overview of the implementation planning is given in the [`nf-core Pangenome Pipeline Planning Phase`](https://docs.google.com/presentation/d/1DzHy_fqs_YH6nMIwxzAPLAaz2CQIlR-k7Y-nQKtn6k8/edit#slide=id.p) presentation.

[![GitHub Actions CI Status](https://github.com/nf-core/pangenome/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/pangenome/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/pangenome/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/pangenome/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/pangenome.svg)](https://hub.docker.com/r/nfcore/pangenome)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23pangenome-4A154B?logo=slack)](https://nfcore.slack.com/channels/pangenome)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/pangenome -profile test,<docker/singularity/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nf-core/pangenome -profile <docker/singularity/conda/institute> --input samplesheet.csv --genome GRCh37
    ```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The nf-core/pangenome pipeline comes with documentation about the pipeline which you can read at [https://nf-core/pangenome/docs](https://nf-core/pangenome/docs) or find in the [`docs/` directory](docs).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits

nf-core/pangenome was originally adapted from the pangenome graph builder [`pggb`](https://github.com/pangenome/pggb) pipeline by Simon Heumos, Michael Heuer.

Many thanks to all who have helped out and contributed along the way, including (but not limited to)\*:

| Name                                                      | Affiliation                                                                           |
|-----------------------------------------------------------|---------------------------------------------------------------------------------------|
| [Philipp Ehmele](https://github.com/imipenem)             | [University of Hamburg, Hamburg, Germany](https://www.uni-hamburg.de/en.html)         |                                         
| [Erik Garisson](https://github.com/ekg)                   | [3Genomics Institute, University of California, Santa Cruz, Santa Cruz, CA, USA, 4Biomolecular Engineering and Bioinformatics, University of California Santa Cruz, Santa Cruz, CA, USA](http://hypervolu.me/~erik/erik_garrison.html)             |                                                                                 
| [Andrea Guarracino](https://github.com/AndreaGuarracino ) | [2University of Rome Tor Vergata, Rome, Italy](http://www.scienze.uniroma2.it/)       |
| [Michael Heuer](https://github.com/heuermh)               | [UC Berkeley, USA](https://https://rise.cs.berkeley.edu)                              |
| [Simon Heumos](https://github.com/subwaystation)          | [QBiC, University of Tübingen, Germany](https://portal.qbic.uni-tuebingen.de/portal/) |
| [Lukas Heumos](https://github.com/zethson)                | [QBiC, University of Tübingen, Germany](https://portal.qbic.uni-tuebingen.de/portal/) |

> \* Listed in alphabetical order


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#pangenome` channel](https://nfcore.slack.com/channels/pangenome) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/pangenome for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
