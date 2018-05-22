'''Snakefile for downloading microarray strand datasets'''
# snakemake -s PlatformStrandFiles.smk

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

RWD = os.getcwd()
HTTP = HTTPRemoteProvider()

with open('scripts/INFORMATION_LISTS/b37_strandfiles.txt', "r") as f:
    platforms = [x.rstrip() for x in f]

rule all:
    input:
        expand('data/platforms/{platforms}/platform.strand', platforms=platforms)


rule StrandFiles:
    input:
        temp(HTTP.remote("http://www.well.ox.ac.uk/~wrayner/strand/{platforms}-strand.zip", keep_local=True))
    output:
        'data/platforms/{platforms}/platform.strand'
    shell:
        "unzip {input} -d data/platforms/{wildcards.platforms} *.strand; "
        "mv data/platforms/{wildcards.platforms}/*.strand data/platforms/{wildcards.platforms}/platform.strand"

#rule StrandFilesMT:
#    input:
#        script
#        path = RWD
#    output:
#        'data/platforms/{platforms}/platform.strand'
#    shell:
#        "unzip {input} -d data/platforms/{wildcards.platforms} *.strand; "
#        "mv data/platforms/{wildcards.platforms}/*.strand data/platforms/{wildcards.platforms}/platform.strand"
