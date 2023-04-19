# ![nf-core/pangenome](docs/images/nf-core-pangenome_logo_light.png#gh-light-mode-only) ![nf-core/pangenome](docs/images/nf-core-pangenome_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/pangenome/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/pangenome)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23pangenome-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/pangenome)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/pangenome** is a bioinformatics best-practice analysis pipeline for pangenome graph construction. The pipeline renders a collection of sequences into a pangenome graph. Its goal is to build a graph that is locally directed and acyclic while preserving large-scale variation. Maintaining local linearity is important for interpretation, visualization, mapping, comparative genomics, and reuse of pangenome graphs.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/pangenome/results).

<p align="center">
    <img title="Pangenome Workflow" src="docs/images/pangenome_workflow.png" width=100%>
</p>

## Pipeline summary

By default, the pipeline currently performs the following:

- All versus all alignment (`WFMASH`)
- Graph induction (`SEQWISH`)
- Graph normalization (`SMOOTHXG`)
- Remove redundancy (`GFAFFIX`)
- Graph statistics and qualitative visualizations (`ODGI`)
- Combine diagnostic information into a report (`MULTIQC`)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run nf-core/pangenome -r dev -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run nf-core/pangenome -r dev --input <BGZIPPED_FASTA> --n_haplotypes <NUM_HAPS_IN_FASTA> --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Advantages over [`PGGB`](https://github.com/pangenome/pggb)

This Nextflow pipeline version's major advantage is that it can distribute the usually computationally heavy all versus all alignment step across a whole cluster. It is capable of splitting the initial approximate alignments into problems of equal size. The base-level alignments are then distributed across several processes. Assuming you have a cluster with 10 nodes and you are the only one using it, we would recommend to set `--wfmash_chunks 10`.
If you have a cluster with 20 nodes, but you have to share it with others, maybe setting it to `--wfmash_chunks 10` could be a good fit, because then you don't have to wait too long for your jobs to finish.

## Documentation

The nf-core/pangenome pipeline comes with documentation about the pipeline [usage](https://nf-co.re/pangenome/usage), [parameters](https://nf-co.re/pangenome/parameters) and [output](https://nf-co.re/pangenome/output).

## Credits

nf-core/pangenome was originally adapted from [PGGB](https://github.com/pangenome/pggb) by [Simon Heumos](https://github.com/subwaystation), [Michael Heuer](https://github.com/heuermh).

> [Simon Heumos](https://github.com/subwaystation) is currently the sole developer.

Many thanks to all who have helped out and contributed along the way, including (but not limited to)\*:

| Name                                                       | Affiliation                                                                                                                                                                                                                                                                                                                                                                       |
| ---------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [Philipp Ehmele](https://github.com/imipenem)              | [Institute of Computational Biology, Helmholtz Zentrum München, Munich, Germany](https://www.helmholtz-muenchen.de/icb/index.html)                                                                                                                                                                                                                                                |
| [Gisela Gabernet](https://github.com/ggabernet)            | [Quantitative Biology Center (QBiC) Tübingen, University of Tübingen, Germany](https://uni-tuebingen.de/en/research/research-infrastructure/quantitative-biology-center-qbic/)                                                                                                                                                                                                    |
| [Erik Garrison](https://github.com/ekg)                    | [The University of Tennessee Health Science Center, Memphis, Tennessee, TN, USA](https://uthsc.edu/)                                                                                                                                                                                                                                                                              |
| [Andrea Guarracino](https://github.com/AndreaGuarracino)   | [Genomics Research Centre, Human Technopole, Milan, Italy](https://humantechnopole.it/en/)                                                                                                                                                                                                                                                                                        |
| [Friederike Hanssen](https://github.com/FriederikeHanssen) | [Quantitative Biology Center (QBiC) Tübingen, University of Tübingen, Germany](https://uni-tuebingen.de/en/research/research-infrastructure/quantitative-biology-center-qbic/) <br> [Biomedical Data Science, Department of Computer Science, University of Tübingen, Germany](https://uni-tuebingen.de/en/faculties/faculty-of-science/departments/computer-science/department/) |
| [Michael Heuer](https://github.com/heuermh)                | [UC Berkeley, USA](https://rise.cs.berkeley.edu)                                                                                                                                                                                                                                                                                                                                  |
| [Lukas Heumos](https://github.com/zethson)                 | [Institute of Computational Biology, Helmholtz Zentrum München, Munich, Germany](https://www.helmholtz-muenchen.de/icb/index.html) <br> [Institute of Lung Biology and Disease and Comprehensive Pneumology Center, Helmholtz Zentrum München, Munich, Germany](https://www.helmholtz-muenchen.de/ilbd/the-institute/cpc/index.html)                                              |
| [Simon Heumos](https://github.com/subwaystation)           | [Quantitative Biology Center (QBiC) Tübingen, University of Tübingen, Germany](https://uni-tuebingen.de/en/research/research-infrastructure/quantitative-biology-center-qbic/) <br> [Biomedical Data Science, Department of Computer Science, University of Tübingen, Germany](https://uni-tuebingen.de/en/faculties/faculty-of-science/departments/computer-science/department/) |
| [Susanne Jodoin](https://github.com/SusiJo)                | [Quantitative Biology Center (QBiC) Tübingen, University of Tübingen, Germany](https://uni-tuebingen.de/en/research/research-infrastructure/quantitative-biology-center-qbic/)                                                                                                                                                                                                    |
| [Júlia Mir Petrol](https://github.com/mirpedrol)           | [Quantitative Biology Center (QBiC) Tübingen, University of Tübingen, Germany](https://uni-tuebingen.de/en/research/research-infrastructure/quantitative-biology-center-qbic/)                                                                                                                                                                                                    |

> \* Listed in alphabetical order

## Acknowledgments

- [QBiC](https://www.qbic.uni-tuebingen.de)
- [deNBI](https://www.denbi.de/)
- [Human Pangenome Reference Consortium](https://humanpangenome.org)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#pangenome` channel](https://nfcore.slack.com/channels/pangenome) (you can join with [this invite](https://nf-co.re/join/slack)), or contact me [Simon Heumos](mailto:simon.heumos@qbic.uni-tuebingen.de?subject=[GitHub]%20nf-core/pangenome).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/pangenome for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## Changelog

[CHANGELOG](CHANGELOG.md)
