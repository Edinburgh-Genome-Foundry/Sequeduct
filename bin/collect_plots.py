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

# Collect plots for easy overview:
import os
import shutil

collected_plot_path = "collected_unfiltered_nanoplots"
nanoplot_path = "results/dir1_preview"
plot_filename = "LengthvsQualityScatterPlot_kde.png"

for directory in os.listdir(nanoplot_path):
    source = os.path.join(nanoplot_path, directory, plot_filename)
    destination = os.path.join(collected_plot_path, directory + "_" + plot_filename)
    shutil.copyfile(source, destination)
