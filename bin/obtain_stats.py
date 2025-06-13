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

preview_nanostat = sys.argv[1]  # skip first arg filename
barcode = sys.argv[2]

na_msg = "FILE NOT FOUND,,"  # commas for CSV
error_msg = "ERROR"  # unexpected structure in NanoStats.txt

def get_stats(filepath):
    try:
        with open(filepath, 'r') as f:
            read_data = f.read().splitlines()
    except:
        return na_msg

    # Number of reads
    number_line = read_data[5]  # should be on this line, but will check below
    if "Number of reads" in number_line:
        # value is at the end of line, remove formatting:
        number_str = number_line.split(" ")[-1].replace(",", "")  # comma is thousand separator
    else:
        number_str = error_msg  # useful for debugging

    # Median read length
    median_line = read_data[3]  # should be on this line
    if "Median read length" in median_line:
        median_str = median_line.split(" ")[-1].replace(",", "")
    else:
        median_str = error_msg

    full_text = ",".join([barcode, number_str, median_str])

    return full_text

full_text = get_stats(preview_nanostat)

print(full_text)
