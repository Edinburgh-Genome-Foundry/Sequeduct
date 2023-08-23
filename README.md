<p align="center">
<img alt="Sequeduct logo" title="Sequeduct" src="images/logo.png" width="120">
</p>

# Sequeduct

![version](https://img.shields.io/badge/current_version-0.2.3-blue)

Sequencing analysis pipeline (aqueduct) for validating plasmids and DNA assembly constructs, using long reads.

## Usage

### Setup

Install [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/).

Pull the Nextflow pipeline:

```bash
nextflow pull edinburgh-genome-foundry/Sequeduct -r v0.2.3
```

Pull the Docker image that contains the required software (requires access to EGF's container repo):

```bash
docker pull ghcr.io/edinburgh-genome-foundry/sequeduct:latest
```

Alternatively, build the image locally from the cloned repo:

```bash
docker build . -f containers/Dockerfile --tag sequeduct_local
```

### Run

Create a directory for your project and copy (or link) the FASTQ directories from your Nanopore run (e.g. into `fastq`). Specify this together with a sample sheet in your commands:

```bash
# Preview
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.2.3 -entry preview --fastq_dir='fastq' \
    --reference_dir='genbank' \
    --sample_sheet='sample_sheet.csv' \
    -profile docker
# Analysis
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.2.3 -entry analysis --fastq_dir='fastq' \
    --reference_dir='genbank' \
    --sample_sheet='sample_sheet.csv' \
    --projectname='EGF project' \
    -profile docker
# Review
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.2.3 -entry review --reference_dir='genbank' \
    --results_csv='results_sheet.csv' \
    --projectname='EGF project review' \
    --all_parts='parts_fasta/part_sequences.fasta' \
    --assembly_plan='assembly_plan.csv' \
    -profile docker
# De novo assembly
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.2.3 -entry assembly --fastq_dir='fastq' \
    --results_csv='results_sheet.csv' \
    --projectname='EGF assembly' \
    -profile docker 
```

The above commands each output a directory within `results`. Similarly, Nextflow creates and uses a directory named `work`, so ensure that your project directory doesn't have one. Specify revision of the project with `-r` (a git branch or tag), and choose a configuration profile (with `-profile`). Profiles are specified in the Nextflow config files.

Use `-with-docker sequeduct_local` to use a locally built Docker image (instead of `-profile docker`).

### Details

For simplicity, the names in the sample sheet are used for finding the reference Genbank files, therefore sample names must match filenames with a ".gb" extension.

Note that canu v2.2 requires minimum 100 reads, otherwise it returns an error. A [fix has been posted](https://github.com/marbl/canu/issues/2035), but it's not released yet.

For convenience, a script is included to collect plot files from the result directories (`bin/collect_plots.py`).

## License = GPLv3+

Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh

Sequeduct was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp), and is released under the GPLv3 license.
