#!/usr/bin/env python
# Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of Sequeduct.
#
# Sequeduct is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Sequeduct is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Sequeduct. If not, see <https:www.gnu.org/licenses/>.

import decimal
import sys

preview_nanostat = sys.argv[1]  # skip first arg filename
analysis_nanostat = sys.argv[2]

na_msg = "NA"
error_msg = "ERROR"

def get_total_count(filepath):
    try:
        with open(filepath, 'r') as f:
            read_data = f.read().splitlines()
    except:
        return na_msg
    
    total_line = read_data[8]  # should be on this line, but will check below
    if "Total bases" not in total_line:
        return error_msg  # useful for debugging
    # value is at the end of line, remove formatting:
    total_str = total_line.split(" ", )[-1].replace(",", "")
    # read count should be an integer:
    total_int = int(float(total_str))
    
    return total_int

preview_bases = get_total_count(preview_nanostat)
analysis_bases = get_total_count(analysis_nanostat)

if type(preview_bases) == int and type(analysis_bases) == int:
    pct = analysis_bases / preview_bases * 100
    pct_str = str(decimal.Decimal(pct).quantize(decimal.Decimal('1'), rounding=decimal.ROUND_HALF_UP))
    print(pct_str + "%", end="")
elif preview_bases == na_msg:  # preview pipeline was not run
    print(na_msg, end="")
else:  # base count could not be obtained
    print(error_msg, end="")
