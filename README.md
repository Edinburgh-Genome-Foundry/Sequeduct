<p align="center">
<img alt="Sequeduct logo" title="Sequeduct" src="images/logo.png" width="120">
</p>

# Sequeduct

![version](https://img.shields.io/badge/current_version-0.1.1-blue)

Sequencing analysis pipeline (aqueduct) for validating plasmids and DNA assembly constructs, using long reads.

## Usage

### Setup

Install [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/).

Download the pipeline:

```bash
git clone git@github.com:Edinburgh-Genome-Foundry/Sequeduct.git
```

Pull the Docker image that contains the required software:

```bash
docker pull ghcr.io/edinburgh-genome-foundry/sequeduct:latest
```

### Run

Create a directory for your project and copy (or link) the FASTQ directories from your Nanopore run (e.g. into `fastq`). Specify this together with a sample sheet in your commands:

```bash
# Preview
nextflow run main.nf -entry preview --fastq_dir='fastq' --reference_dir='genbank' \
    --sample_sheet='sample_sheet.csv' \
    -with-docker ghcr.io/edinburgh-genome-foundry/sequeduct
# Analysis
nextflow run main.nf -entry analysis --fastq_dir='fastq' --reference_dir='genbank' \
    --sample_sheet='sample_sheet.csv' \
    --projectname='EGF project' \
    -with-docker ghcr.io/edinburgh-genome-foundry/sequeduct
# Review
nextflow run main.nf -entry review --reference_dir='genbank' \
    --results_csv='results_finalised.csv' \
    --projectname='EGF project review' \
    --all_parts='parts_fasta/part_sequences.fasta' \
    --assembly_plan='assembly_plan.csv' \
    -with-docker ghcr.io/edinburgh-genome-foundry/sequeduct
```

The above three commands each output a directory in `results`. Similarly, NextFlow creates and uses a directory named `work`, so ensure that your project directory doesn't have one.

## Build the image locally

```bash
git clone git@github.com:Edinburgh-Genome-Foundry/Ediacara.git
docker build . -f containers/Dockerfile --tag sequeduct_local
```

## Copyright

Copyright 2021 Edinburgh Genome Foundry

Sequeduct was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp).
