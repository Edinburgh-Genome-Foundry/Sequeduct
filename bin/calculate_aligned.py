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

fastq_file = sys.argv[1]  # skip first arg filename
paf_file = sys.argv[2]

import decimal

from Bio import SeqIO

import ediacara

# FILTERED READS
fastq_list = []
with open(fastq_file) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        fastq_list += [record.id]
fastq_readset = set(fastq_list)

# ALIGNED READS
paf = ediacara.ComparatorGroup.load_paf(paf_file)
aligned_readset = set(paf["query_name"])

# DIFFERENCE
pct = len(aligned_readset) / len(fastq_readset)  * 100
pct_str = str(decimal.Decimal(pct).quantize(decimal.Decimal('1'), rounding=decimal.ROUND_HALF_UP))
print(pct_str, end="")
