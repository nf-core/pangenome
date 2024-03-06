# nf-core/pangenome: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Input Check](#input-check)
  - [bgzip](#bgzip)
  - [samtools faidx](#samtools-faidx)
- [Community Detection](#community-detection)
  - [paf2net](#paf2net)
  - [net2communities](#net2communities)
  - [extract communities](#extract-communities)
- [wfmash](#wfmash)
  - [wfmash map community](#wfmash-map-community)
  - [wfmash map](#wfmash-map)
  - [wfmash align](#wfmash-align)
  - [split approx mappings in chunks](#split-approx-mappings-in-chunks)
- [seqwish](#seqwish)
- [smoothxg](#smoothxg)
- [gfaffix](#gfaffix)
- [odgi](#odgi)
  - [odgi build](#odgi-build)
  - [odgi stats](#odgi-stats)
  - [odgi sort](#odgi-sort)
  - [odgi unchop](#odgi-unchop)
  - [odgi view](#odgi-view)
  - [odgi viz](#odgi-viz)
  - [odgi layout](#odgi-layout)
  - [odgi draw](#odgi-draw)
  - [odgi squeeze](#odgi-squeeze)
- [vg](#vg)
  - [vg deconstruct](#vg-deconstruct)
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline -> final report(s)!
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Input Check

### bgzip

[bgzip](http://www.htslib.org/doc/bgzip.html) compresses a file into a series of small 'BGZF' blocks in a similar manner to, and compatible with, gzip. This allows indexes to be built against the compressed file and used to retrieve portions of the data without having to decompress the entire file.

<details markdown="1">
<summary>Output files</summary>

- `tabix_bgzip/`
  - `<INPUT_FASTA>.gz`: The bgzip compressed input FASTA file.
  - `<INPUT_FASTA>.community.[0-9]{1,}.fa.gz`: A bgzipped compressed community FASTA file. _Only appears when `--communities` is provided._
  </details>

### samtools faidx

[samtools faidx](http://www.htslib.org/doc/samtools-faidx.html) indexes or queries regions from a FASTA file.

<details markdown="1">
<summary>Output files</summary>

- `samtools_faidx/`
  - `<INPUT_FASTA>.fai`: FASTA index of the file provided via `--input`.
  - `<INPUT_FASTA>.gzi`: Compressed FASTA index of the file provided via `--input`.
  - `<INPUT_FASTA>.community.[0-9]{1,}.fa.gz.fai`: FASTA index of a community FASTA file. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.community.[0-9]{1,}.fa.gz.gzi`: Compressed FASTA index of a community FASTA file. _Only appears when `--communities` is provided._
  </details>

<!-- [MultiQC - TEST](images/pangenome_workflow.png) -->

## Community Detection

### paf2net

[paf2net](https://github.com/pangenome/pggb/blob/master/scripts/paf2net.py) is python script that projects wfmash's PAF mappings (the implied overlap and containment graph) into an edge list, a list of edge weights, and an 'id to sequence name' map.

<details markdown="1">
<summary>Output files</summary>

- `paf2net/`
  - `<INPUT_FASTA>.paf.vertices.id2name.txt`: TXT file with a mapping of vertex identifiers to sequence names. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.paf.edges.weights.txt`: TXT file with the weights of the edges. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.paf.edges.list.txt`: TXT file listing all edge connections. _Only appears when `--communities` is provided._

</details>

### net2communities

[net2communities](https://github.com/pangenome/pggb/blob/master/scripts/net2communities.py) is a python script that detects communities by applying the Leiden algorithm ([Trag et al., Nature 2019](https://www.nature.com/articles/s41598-019-41695-z)).

<details markdown="1">
<summary>Output files</summary>

- `net2communities/`
  - `<INPUT_FASTA>.community.[0-9](1,).txt`: TXT file with the sequence names of the community. _Only appears when `--communities` is provided._
  </details>

### extract communities

[extract communities](https://github.com/nf-core/pangenome/blob/dev/modules/local/extract_communities/main.nf) is a locally modified module of [samtools faidx](#samtools-faidx). The original module only allowed for the indexing of FASTA files, but not for the extraction of single sequences.

<details markdown="1">
<summary>Output files</summary>

- `extract_communities/`
  - `<INPUT_FASTA>.community.[0-9](1,).txt.fa`: A community FASTA file. _Only appears when `--communities` is provided._
  </details>

## wfmash

[wfmash](https://github.com/waveygang/wfmash) is an aligner for pangenomes based on sparse homology mapping and wavefront inception. In this pipeline it is used in various processes and ways.

<details markdown="1">
<summary>Output files</summary>

- `wfmash/`
  - `<INPUT_FASTA>.paf`: PAF file containing CIGAR strings of all pairwise alignments.
  - `<INPUT_FASTA>.community.[0-9]{1,}.paf`: Community PAF file containing CIGAR strings of all pairwise alignments of the specific community. _Only appears when `--communities` is provided._
  </details>

### wfmash map community

Here `wfmash` was applied in approximate mappping mode in order to create all-all alignments of all input sequences for their subsequent community discovery.

<details markdown="1">
<summary>Output files</summary>

- `wfmash_map_community/`
  - `<INPUT_FASTA>.paf`: PAF file containing approximate mappings of all pairwise alignments.
  </details>

### wfmash map

Here `wfmash` was applied in approximate mappping mode in order to create all-all alignments of all input sequences. We can then split the base pair level alignment problem into several equal problem sizes with [split approx mappings in chunks](#split-approx-mappings-in-chunks).

<details markdown="1">
<summary>Output files</summary>

- `wfmash_map/`
  - `<INPUT_FASTA>.paf`: PAF file containing approximate mappings of all pairwise alignments.
  - `<INPUT_FASTA>.community.[0-9]{1,}.paf`: Community PAF file containing approximate mappings of all pairwise alignments of the specific community. _Only appears when `--communities` is provided._
  </details>

### split approx mappings in chunks

[split approx mappings in chunks](https://github.com/waveygang/wfmash/blob/master/scripts/split_approx_mappings_in_chunks.py) is a python script that takes the approximate mappings, weighs each mapping by computing its length \* (1 - estimated identity), then creates N new files where the mapping sets have a similar sum of weights.

<details markdown="1">
<summary>Output files</summary>

- `wfmash_map/`
  - `<INPUT_FASTA>.paf.chunk_[0-9]{1,}.paf`: PAF file containing base level alignments of a specific chunk.
  - `<INPUT_FASTA>.community.[0-9]{1,}.paf.chunk_[0-9]{1,}.paf`: Community PAF file containing base level alignments of a specific chunk of the specific community. _Only appears when `--communities` is provided._
  </details>

### wfmash align

Here `wfmash` was applied in base pair level alignment mode in order to refine the approximate all-all alignments of all input sequences.

<details markdown="1">
<summary>Output files</summary>

- `wfmash_map/`
  - `<INPUT_FASTA>.chunk_[0-9]{1,}.paf`: PAF file containing base level alignments of a specific chunk.
  - `<INPUT_FASTA>.community.[0-9]{1,}.chunk_[0-9]{1,}.paf`: Community PAF file containing base level alignments of a specific chunk of the specific community. _Only appears when `--communities` is provided._
  </details>

## seqwish

[seqwish](https://github.com/ekg/seqwish) implements a lossless conversion from pairwise alignments between sequences to a variation graph encoding the sequences and their alignments.

<details markdown="1">
<summary>Output files</summary>

- `seqwish/`
  - `<INPUT_FASTA>.seqwish.gfa`: Raw pangenome graph induced from the all-versus-all alignments.
  - `<INPUT_FASTA>.community.[0-9]{1,}.seqwish.gfa`: Community GFA file containing the raw pangenome graph induced from the all-versus-all alignments of the specific community. _Only appears when `--communities` is provided._
  </details>

## smoothxg

[smoothxg](https://github.com/pangenome/smoothxg) finds blocks of paths that are collinear within a variation graph. It applies partial order alignment to each block, yielding an acyclic variation graph. Then, to yield a "smoothed" graph, it walks the original paths to lace these subgraphs together. The resulting graph only contains cyclic or inverting structures larger than the chosen block size, and is otherwise manifold linear.

<details markdown="1">
<summary>Output files</summary>

- `smoothxg/`
  - `<INPUT_FASTA>.smoothxg.gfa`: Smoothed pangenome graph in GFA format.
  - `<INPUT_FASTA>.community.[0-9]{1,}.smoothxg.gfa`: Community GFA file containing the smoothed pangenome graph of the specific community. _Only appears when `--communities` is provided._
  </details>

## gfaffix

[gfaffix](https://github.com/marschall-lab/GFAffix) identifies walk-preserving shared affixes in variation graphs and collapses them into a non-redundant graph structure.

<details markdown="1">
<summary>Output files</summary>

- `gfaffix/`
  - `<INPUT_FASTA>.gfaffix.gfa`: Non-node-redundant pangenome graph in GFA format.
  - `<INPUT_FASTA>.gfaffix.txt`: Graph shared affixes in TSV format.
  - `<INPUT_FASTA>.community.[0-9]{1,}.gfaffix.gfa`: Community GFA file containing the non-node-redundant pangenome graph of the specific community. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.community.[0-9]{1,}.gfaffix.txt`: Community TSV file containing the graph shared affixes in TSV format of the specific community. _Only appears when `--communities` is provided._
  </details>

## odgi

[odgi](https://github.com/pangenome/odgi) provides an efficient and succinct dynamic DNA sequence graph model, as well as a host of algorithms that allow the use of such graphs in bioinformatic analyses. In this pipeline, a huge variety of odgi's subcommands are used to process the built graphs. As a rule of thumb, all files in ODGI format end with `.og`.

### odgi build

[odgi build](https://odgi.readthedocs.io/en/latest/rst/commands/odgi_build.html) constructs a dynamic succinct variation graph in ODGI format from a GFAv1.

<details markdown="1">
<summary>Output files</summary>

- `odgi_build/`
  - `<INPUT_FASTA>.*.og`: Graph in ODGI format.
  - `<INPUT_FASTA>.community.[0-9]{1,}.*.og`: Community ODGI file of the specific community. _Only appears when `--communities` is provided._
  </details>

### odgi stats

[odgi stats](https://odgi.readthedocs.io/en/latest/rst/commands/odgi_stats.html) describes various metrics of a variation graph.

<details markdown="1">
<summary>Output files</summary>

- `odgi_stats/`
  - `<INPUT_FASTA>.*.og.stats.yaml`: YAML file with graph metrics.
  - `<INPUT_FASTA>.community.[0-9]{1,}.*.og.stats.yaml`: Community YAML file with graph metrics of the specific community. _Only appears when `--communities` is provided._
  </details>

### odgi sort

[odgi sort](https://pangenome.github.io/odgi.github.io/rst/commands/odgi_sort.html) sorts a succinct variation graph. it offers a diverse palette of sorting algorithms to determine the node order.

> In the pipeline itself, this is the last tool invoked before the final graph in ODGI format is written on disk.

<details markdown="1">
<summary>Output files</summary>

- `FINAL_ODGI/`
  - `<INPUT_FASTA>.*.Ygs.og`: Sorted variation graph in ODGI format.
  - `<INPUT_FASTA>.community.[0-9]{1,}.*.Ygs.og`: Community ODGI file with a sorted variation graph of the specific community. _Only appears when `--communities` is provided._
  </details>

The order of the sortings:

- `Y`: PG-SGD
- `g`: grooming
- `s`: topolocigal sort

### odgi unchop

[odgi unchop](https://pangenome.github.io/odgi.github.io/rst/commands/odgi_unchop.html) merges unitigs into single nodes.

<details markdown="1">
<summary>Output files</summary>

- `odgi_unchop/`
  - `<INPUT_FASTA>.*.unchop.og`: Unchopped variation graph in ODGI format.
  - `<INPUT_FASTA>.community.[0-9]{1,}.*.unchop.og`: Community ODGI file with a unchopped variation graph of the specific community. _Only appears when `--communities` is provided._
  </details>

### odgi view

[odgi view](https://pangenome.github.io/odgi.github.io/rst/commands/odgi_view.html) can convert a graph in ODGI format to GFAv1.

<details markdown="1">
<summary>Output files</summary>

- `FINAL_GFA/`
  - `<INPUT_FASTA>.*.view.og`: Final graph in GFAv1 format.
  - `<INPUT_FASTA>.community.[0-9]{1,}.*.view.og`: Final community GFAv1 file the specific community. _Only appears when `--communities` is provided._
  </details>

### odgi viz

[odgi viz](https://odgi.readthedocs.io/en/latest/rst/commands/odgi_viz.html) visualizes a variation graph in 1D.

<details markdown="1">
<summary>Output files</summary>

- `odgi_viz/`
  - `<INPUT_FASTA>.gfaffix.viz_*_multiqc.png`: 1D visualization of a genome variation graph in PNG format ready to be put into a [MultiQC](#multiqc) report.
  - `<INPUT_FASTA>.squeeze.viz_*_multiqc.png`: 1D visualization of all communities combined in one graph in PNG format ready to be put into a [MultiQC](#multiqc) report. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.community.[0-9]{1,}.gfaffix.viz_*_multiqc.png`: 1D visualizaton of a community genome variation graph ready to be put into a [MultiQC](#multiqc) report. _Only appears when `--communities` is provided._
  </details>

### odgi layout

[odgi layout](https://pangenome.github.io/odgi.github.io/rst/commands/odgi_layout.html) uses the PG-SGD algorithm to calculate a 2D layout of a variation graph. The layout in TSV format and a corresponding GFAv1 graph from folder `FINAL_GFA` can be loaded into [waragraph](https://github.com/chfi/waragraph) for an interactive visualization.

<details markdown="1">
<summary>Output files</summary>

- `odgi_layout/`
  - `<INPUT_FASTA>.gfaffix.tsv`: 2D layout in TSV format.
  - `<INPUT_FASTA>.gfaffix.lay`: 2D layout in binary LAY format.
  - `<INPUT_FASTA>.squeeze.tsv`: 2D layout in TSV format of all communities combined in one graph. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.squeeze.tsv`: 2D layout in binary LAY format of all communities combined in one graph. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.community.[0-9]{1,}.gfaffix.tsv`: 2D layout in TSV format of a community genome variation graph. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.community.[0-9]{1,}.gfaffix.tsv`: 2D layout in binary LAY format of a community genome variation graph. _Only appears when `--communities` is provided._
  </details>

### odgi draw

[odgi draw](https://pangenome.github.io/odgi.github.io/rst/commands/odgi_draw.html) takes a 2D graph layout in binary LAY format and a corresponding variation graph in GFAv1 or ODGI format and renders a static 2D visualization.

<details markdown="1">
<summary>Output files</summary>

- `odgi_draw/`
  - `<INPUT_FASTA>.gfaffix.draw_multiqc.png`: 2D visualization of a genome variation graph in PNG format ready to be put into a [MultiQC](#multiqc) report.
  - `<INPUT_FASTA>.gfaffix.png`: Low resolution 2D visualization of a genome variation graph in PNG format.
  - `<INPUT_FASTA>.squeeze.draw_multiqc.png`: 2D visualization of all communities combined in one graph in PNG format ready to be put into a [MultiQC](#multiqc) report. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.squeeze.png`: Low resolution 2D visualization of all communities combined in one graph in PNG format. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.community.[0-9]{1,}.gfaffix.draw_multiqc.png`: 2D visualizaton of a community genome variation graph ready to be put into a [MultiQC](#multiqc) report. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.community.[0-9]{1,}.gfaffix..png`: Low resolution 2D visualizaton of a community genome variation graph. _Only appears when `--communities` is provided._
  </details>

### odgi squeeze

[odgi squeeze](https://pangenome.github.io/odgi.github.io/rst/commands/odgi_squeeze.html) puts multiple variation graphs in one file.

<details markdown="1">
<summary>Output files</summary>

- `FINAL_ODGI/`
  - `<INPUT_FASTA>.squeeze.og`: All graphs of all communities combined in one graph in ODGI format. _Only appears when `--communities` is provided._
  </details>

## vg

[vg](https://github.com/vgteam/vg) is the `v`ariation `g`raph toolkit for data structures, interchange formats, alignment, genotyping, and variant calling methods of genome variation graphs.

### vg deconstruct

[vg deconstruct](https://github.com/vgteam/vg) outputs VCF records for snarls present in a graph relative to one or several chosen reference path(s).

<details markdown="1">
<summary>Output files</summary>

- `vg_deconstruct/`
  - `<INPUT_FASTA>.gfafix.*.vcf`: Variants in VCF format of the graph.
  - `<INPUT_FASTA>.gfafix.*.vcf.stats`: Statistics of the variants in VCF format of the graph.
  - `<INPUT_FASTA>.gfafix.*.decomposed.vcf`: Decomposed variants in VCF format of the graph.
  - `<INPUT_FASTA>.gfafix.*.decomposed.vcf.stats`: Statistics of the decomposed variants in VCF format of the graph.
  - `<INPUT_FASTA>.squeeze.*.vcf`: Variants in VCF format of the graph containing all communities. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.squeeze.*.vcf.stats`: Statistics of the variants in VCF format of the graph containing all communities. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.squeeze.*.decomposed.vcf`: Decomposed variants in VCF format of the graph containing all communities. _Only appears when `--communities` is provided._
  - `<INPUT_FASTA>.squeeze.*.decomposed.vcf.stats`: Statistics of the decomposed variants in VCF format of the graph containing all communities. _Only appears when `--communities` is provided._
  </details>

## MultiQC

In the ODGI table section of the MultiQC report, it can happen that one observes the actual sample name and `seqwish`.
The `seqwish` sample is the graph which was produced by [seqwish](#seqwish).
The named sample, which just contains the sample name in the name is the final graph:

- [seqwish](#seqwish) -> [smoothxg](#smoothxg) -> [gfaffix](#gfaffix) -> [odgi build](odgi-build) -> [odgi unchop](#odgi-unchop) -> [odgi sort](#odgi-sort)

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
