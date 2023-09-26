<p align="center">
<img alt="Sequeduct logo" title="Sequeduct" src="images/logo.png" width="120">
</p>

# Sequeduct

![version](https://img.shields.io/badge/current_version-0.3.1-blue)

Sequencing analysis pipeline (aqueduct) for validating plasmids and DNA assembly constructs, using long reads.

## Usage

### Setup

Install [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/).

Pull the Nextflow pipeline:

```bash
nextflow pull edinburgh-genome-foundry/Sequeduct -r v0.3.1
```

#### Docker image

Build the image that contains the software required for running the pipeline. First, obtain the code (Dockerfile) either by downloading or cloning:

##### Download

Download the repository...

* click on the "<> Code" button at the top of this page, and 'Download ZIP'
* open a terminal where the file was downloaded
* Unzip the file (e.g. `unzip Sequeduct-main.zip`)

##### Clone

... or clone the repository:

```bash
git clone https://github.com/Edinburgh-Genome-Foundry/Sequeduct.git
```

#### Build

Change to the downloaded directory (e.g. `cd Sequeduct-main/`), then run:

```bash
docker build . -f containers/Dockerfile --tag sequeduct_local
```

where sequeduct_local is a custom tag that you can specify, and should be used in the run commands below.

Alternatively, pull the Docker image if you have access to EGF's container repo (e.g. EGF staff members):

```bash
docker pull ghcr.io/edinburgh-genome-foundry/sequeduct:v0.3.1
```

Use `-profile docker` to use this image in the below commands, instead of `-with-docker sequeduct_local`.


### Run

Create a directory for your project and copy (or link) the FASTQ directories from your Nanopore run (e.g. `fastq_pass`). Specify this together with a sample sheet in your commands:

```bash
# Preview
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.3.1 -entry preview --fastq_dir='fastq_pass' \
    --reference_dir='genbank' \
    --sample_sheet='sample_sheet.csv' \
    -with-docker sequeduct_local
# Analysis
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.3.1 -entry analysis --fastq_dir='fastq_pass' \
    --reference_dir='genbank' \
    --sample_sheet='sample_sheet.csv' \
    --projectname='EGF project' \
    -with-docker sequeduct_local
# Review
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.3.1 -entry review --reference_dir='genbank' \
    --results_csv='results_sheet.csv' \
    --projectname='EGF project review' \
    --all_parts='parts_fasta/part_sequences.fasta' \
    --assembly_plan='assembly_plan.csv' \
    -with-docker sequeduct_local
# De novo assembly
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.3.1 -entry assembly --fastq_dir='fastq_pass' \
    --assembly_sheet='assembly_sheet.csv' \
    -with-docker sequeduct_local
```

The above commands each output a directory within a created `results` directory. Similarly, Nextflow creates and uses a directory named `work`, so ensure that your project directory doesn't have a directory with the same name. Specify revision of the project with `-r` (a git branch or tag), and choose a configuration profile (with `-profile`). Profiles are specified in the Nextflow config files. The Review pipeline utilises the output files of the Analysis pipeline, but otherwise the pipelines are independent. Please find example sheets in the `examples` directory.

A more detailed example and demonstration data are available at the [Sequeduct demo](https://github.com/Edinburgh-Genome-Foundry/Sequeduct_demo) site.


### Details

For simplicity, the names in the sample sheet are used for finding the reference Genbank files, therefore sample names must match filenames with a ".gb" extension.

If you have the FASTQ files in gzip compressed format (`.gz`), then you must uncompress them (e.g. run `gunzip --recursive *` in the FASTQ folder).

Note that canu v2.2 requires minimum 100 reads, otherwise it returns an error. A [fix has been posted](https://github.com/marbl/canu/issues/2035), but it's not released yet.

For convenience, a script is included to collect plot files from the result directories (`bin/collect_plots.py`).

The pipeline was designed to work with data from one or more barcodes (FASTQ subdirectories). It has been tested on a desktop machine running Ubuntu 20.04.6 LTS (Memory: 15.5 GiB; CPU: Intel® Core™ i5-6500 CPU @ 3.20GHz × 4), and confirmed to work with up to 96 barcodes. The largest tested dataset was 1.5 GB Nanopore FASTQ data, resulting in 1.1 GB filtered data (100k filtered reads) with up to 55 MB individual filtered FASTQ files (i.e. per sample). If the dataset is much larger, then it may return an error at the variant call or another step. A recommended solution is to increase the quality cutoff (with parameter `--quality_cutoff`), and optionally the minimum length cutoff (`--min_length`), to work with fewer but better reads.


## License = GPLv3+

Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh

Sequeduct was designed by [Giovanni Stracquadanio](https://github.com/stracquadaniolab/) and Peter Vegh. It's implemented in Nextflow by [Peter Vegh](https://github.com/veghp) at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/), and is released under the GPLv3 license.
