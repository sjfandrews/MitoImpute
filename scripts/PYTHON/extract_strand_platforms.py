#!/usr/bin/env python

SF="/Users/TimMcInerney/GitCode/MitoImpute/scripts/INFORMATION_LISTS/b37_strandfiles.txt"
outf="/Users/TimMcInerney/GitCode/MitoImpute/scripts/INFORMATION_LISTS/b37_platforms.txt"

with open(SF, "r") as f:
    platforms = [x.rstrip() for x in f]

with open(outf, 'wr') as o:
    for i in platforms:
        o.write(i + "\n")