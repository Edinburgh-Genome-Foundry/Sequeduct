#!/usr/bin/env python
import sys

genbank_path = sys.argv[1]  # skip first filename
sample = sys.argv[2]
sample_fasta = sys.argv[3]
max_len_fraction = sys.argv[4]
flag = sys.argv[5]

from Bio import SeqIO

# Genbank in
record = SeqIO.read(genbank_path, "genbank")
record.id = sample
# FASTA out
with open(sample_fasta, "w") as output_handle:
    SeqIO.write(record, output_handle, "fasta")

if flag == "none":
    pass  # no need for stdout
elif flag == "length":
    cutoff_length = int(float(len(record)) * float(max_len_fraction))
    print(cutoff_length, end="")
elif flag == "canu":
    print(str(round(len(record) / 1000)), end="")  # to get k value for canu
