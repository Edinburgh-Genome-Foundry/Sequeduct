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

import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

fasta_filename = sys.argv[1]  # skip first (filename)
genbank_files = sys.argv[2:]

open_fasta = open(fasta_filename, "a")

for genbank in genbank_files:
    with open(genbank, "r") as handle:
        record = SeqIO.read(
            handle,
            "genbank",
        )
        name = Path(genbank).stem
        sequence = str(record.seq)
        record_2x = SeqRecord(Seq(sequence + sequence))
        record_2x.id = name
        record_2x.name = name
        record_2x.description = ""

        SeqIO.write(record_2x, open_fasta, "fasta")

open_fasta.close()
