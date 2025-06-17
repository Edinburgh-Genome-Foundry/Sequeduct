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

import pandas as pd

preview_csv = sys.argv[1]  # skip first arg filename
sorted_csv = sys.argv[2]

preview_df = pd.read_csv(preview_csv, header=None)
preview_df.columns = ["Barcode", "Number_of_reads", "Median_read_length"]
preview_df.sort_values(by=["Barcode"], inplace=True)

print(sorted_csv)

preview_df.to_csv(sorted_csv, index=False)
