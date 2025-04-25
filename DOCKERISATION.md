# Docker image

Install [Docker](https://www.docker.com/).

Build the image that contains the software required for running the pipeline. First, obtain the code (Dockerfile) either by downloading or cloning:

### Download

Download the repository...

* click on the "<> Code" button at the top of this page, and 'Download ZIP'
* open a terminal where the file was downloaded
* Unzip the file (e.g. `unzip Sequeduct-main.zip`)

### Clone

... or clone the repository:

```bash
git clone https://github.com/Edinburgh-Genome-Foundry/Sequeduct.git
```

## Build

Change to the downloaded directory (e.g. `cd Sequeduct-main/`), then run:

```bash
docker build . -f containers/Dockerfile --tag sequeduct_local
```

where `sequeduct_local` is a custom tag that you can specify, and should be used in the run commands. For example:

```bash
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.4.3 -entry analysis --fastq_dir='fastq_pass' \
    --reference_dir='genbank' \
    --sample_sheet='sample_sheet.csv' \
    --projectname='EGF project' \
    -with-docker sequeduct_local
```

## With access to EGF's container repo

Alternatively, pull the Docker image if you have access to EGF's container repo (e.g. EGF staff members):

```bash
docker pull ghcr.io/edinburgh-genome-foundry/sequeduct:v0.4.3
```

Use `-profile docker` to use this image. Example:

```bash
nextflow run edinburgh-genome-foundry/Sequeduct -r v0.4.3 -entry analysis --fastq_dir='fastq_pass' \
    --reference_dir='genbank' \
    --sample_sheet='sample_sheet.csv' \
    --projectname='EGF project' \
    -profile docker
```
