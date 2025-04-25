#!/usr/bin/env python
# Copyright 2025 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of Sequeduct.
#
# Sequeduct is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Sequeduct is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Sequeduct. If not, see <https:www.gnu.org/licenses/>.
import os
import sys
from Bio import SeqIO
import ediacara as edi

split_fastq_dir = sys.argv[1]  # skip first (filename)
barcode = sys.argv[2]
paf_file = sys.argv[3]
fastq_files = sys.argv[4:]

paf = edi.ComparatorGroup.load_paf(paf_file)

# dict = {ref : [reads] , ... }
ref_read_dict = {target: [] for target in paf.target_name.unique()}
for index, row in paf.iterrows():
    ref_read_dict[row["target_name"]] += [row["query_name"]]


def invert_dict(d):
    inverse_dict = dict()
    for ref, reads in d.items():
        for read in reads:
            if read not in inverse_dict:
                inverse_dict[read] = [ref]
            else:
                inverse_dict[read].append(ref)
    return inverse_dict


read_ref_dict = invert_dict(ref_read_dict)

# Remove reads with alignments to more than one ref:
read_ref_filtered_dict = dict()
for read, refs in read_ref_dict.items():
    dedup_ref = list(set(refs))
    if len(dedup_ref) == 1:
        read_ref_filtered_dict[read] = dedup_ref[0]

# Split fastq using Biopython:
ref_dir_dict = dict()
counter = 1
for ref in paf.target_name.unique():
    dirname = os.path.join(split_fastq_dir, barcode + "_" + str(counter))
    ref_dir_dict[ref] = dirname
    try:
        os.mkdir(dirname)
    except:
        pass
    counter += 1

ref_file_dict = {
    ref: os.path.join(directory, ref + ".fastq")
    for ref, directory in ref_dir_dict.items()
}

print(ref_file_dict)

os.mkdir(split_fastq_dir)
for subdir in ref_dir_dict.values():
    os.mkdir(subdir)

for ref, out_fastq in ref_file_dict.items():
    ref_file_dict[ref] = open(out_fastq, "a")

for in_fastq_file in fastq_files:
    with open(in_fastq_file, "r") as handle:
        for fastq_read in SeqIO.parse(handle, "fastq"):
            try:
                ref = read_ref_filtered_dict[fastq_read.id]
                SeqIO.write(fastq_read, ref_file_dict[ref], "fastq")
            except:  # read was filtered out
                pass

for filehandle in ref_file_dict.values():
    filehandle.close()
