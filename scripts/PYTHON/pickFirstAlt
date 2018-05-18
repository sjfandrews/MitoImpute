#!/usr/bin/env python3
from pysam import VariantFile
import sys

vcf_in = VariantFile(sys.argv[1], 'r')
vcf_out = VariantFile('-', 'w', header=vcf_in.header)
cp = (0, 0)
for rec in vcf_in.fetch():
    if (rec.chrom, rec.pos) != cp:
        vcf_out.write(rec)
    cp = (rec.chrom, rec.pos)
