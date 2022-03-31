#!/usr/bin/env python
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
