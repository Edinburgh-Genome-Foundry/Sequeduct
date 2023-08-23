#!/usr/bin/env python
# Copyright 2021 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of Sequeduct.
#
# Sequeduct is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Sequeduct is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Sequeduct. If not, see <https:www.gnu.org/licenses/>.

import sys

assembly_dir = sys.argv[1]  # skip first filename
params_assembly_prefix = sys.argv[2]
params_canu_postfix = sys.argv[3]
trimmed_denovo = sys.argv[4]
barcode = sys.argv[5]

from Bio import SeqIO

canu_fasta = assembly_dir + '/' + params_assembly_prefix + params_canu_postfix
try:
    contig = SeqIO.read(canu_fasta, format="fasta")
except:
    print("The FASTA file contains more than 1 contigs. First contig used.")
    contig = next(SeqIO.parse(canu_fasta, format="fasta"))

entries = contig.description.split(" ")
desc_dict = {"name": entries[0]}  # first is the name
for entry in entries[1:]:  # addressed the first one above
    elements = entry.split("=")
    desc_dict[elements[0]] = elements[1]

# canu assembly: 0-based, from-index inclusive, end-index exclusive
if desc_dict["suggestCircular"] == "yes":  # as output by canu
    start, end = desc_dict["trim"].split("-")  # must contain 2 values
    start = int(start)
    end = int(end)
    SeqIO.write(contig[start:end], trimmed_denovo, format="fasta")
else:  # keep intact
    SeqIO.write(contig, trimmed_denovo, format="fasta")

print("Trimmed:", barcode)
