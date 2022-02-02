<p align="center">
<img alt="Sequeduct logo" title="Sequeduct" src="images/logo.png" width="120">
</p>

# Sequeduct

Sequencing analysis pipeline (aqueduct) for validating plasmids and DNA assembly constructs, using long reads.

## Build

```bash
cd containers
docker build --tag sequeduct .
```

## Usage

```bash
# Preview
nextflow run sequeduct_preview.nf --fastq_dir='fastq' --sample_sheet='sample_sheet.csv'
# Analysis
nextflow run sequeduct_analysis.nf --fastq_dir='fastq' --reference_dir='ref' --sample_sheet='sample_sheet.csv' --projectname='EGF project' -with-docker sequeduct
# Review
nextflow run /home/peter/data/lib/software/sequeduct_software/Sequeduct/sequeduct_review.nf --reference_dir='ref' --results_csv='results.csv' --projectname='EGF project' --all_parts='part_sequences.fasta' --assembly_plan='assembly_plan.csv'
```

## Copyright

Copyright 2021 Edinburgh Genome Foundry

Sequeduct was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp).
