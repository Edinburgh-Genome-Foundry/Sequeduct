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
import re
import cyvcf2

vcf_path = sys.argv[1]  # skip first filename
out_vcf = sys.argv[2]

vcf = cyvcf2.VCF(vcf_path)
w = cyvcf2.Writer(out_vcf, vcf)

for variant in vcf:
    # cf. Ediacara
    skip = False
    result = re.search(r"((\w)\2{4,})", variant.REF)
    if result is not None:
        skip = True
    else:  # no need to check this if regex was positive
        # ignore if ALT < 50%: mutation is not real or sample is polymorphic
        if variant.INFO.get("AO") < variant.INFO.get("DP") * 0.5:
            skip = True

    if not skip:
        w.write_record(variant)

w.close()
vcf.close()
