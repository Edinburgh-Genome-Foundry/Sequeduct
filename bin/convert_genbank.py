#!/usr/bin/env python
import sys

genbank_path = sys.argv[1]  # skip first filename
sample = sys.argv[2]
sample_fasta = sys.argv[3]
flag = sys.argv[4]

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
    print(len(record), end="")
elif flag == "canu":
    print(str(round(len(record) / 1000)), end="")  # to get k value for canu
