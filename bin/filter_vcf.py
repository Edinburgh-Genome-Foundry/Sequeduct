#!/usr/bin/env python
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
