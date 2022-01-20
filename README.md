<p align="center">
<img alt="Sequeduct logo" title="Sequeduct" src="images/logo.png" width="120">
</p>


# Sequeduct

Sequencing analysis pipeline (aqueduct) for validating plasmids and DNA assembly constructs, using long reads.


## Usage

```bash
nextflow run sequeduct_preview.nf --fastq_dir='fastq' --sample_sheet='sample_sheet.csv'

nextflow run sequeduct_analysis.nf --fastq_dir='fastq' --references='ref' --sample_sheet='sample_sheet.csv'
```


## Copyright

Copyright 2021 Edinburgh Genome Foundry

Sequeduct was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp).
