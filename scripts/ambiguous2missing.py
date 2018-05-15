#!/usr/bin/env python2
## IMPORT NECESSARY PACKAGES
import os
import argparse
import numpy
from tqdm import *
import time

## DEFINE MAIN SCRIPT/FUNCTION
def main():
    #START THE TIMER
    start_time = time.time()
    #DEFINE INPUT ARGUMENTS
    parser=argparse.ArgumentParser(description='ambiguous2missing',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile', type=str, required=True, help='input fasta file (.fasta)')
    parser.add_argument('-o', '--outfile', dest='outfile', type=str, required=False, help='output fasta file (.fasta)')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", required=False, help='turn on verbose mode')

    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile
    verbose = args.verbose

    if outfile is None:
        outfile = os.path.splitext(infile)[0] + '_ambig2missing.fasta'

    if verbose:
        print('\nFASTA FILE IN:\t\t{}'.format(infile))
        print('FASTA FILE IN:\t\t{}\n'.format(outfile))

    #DISPLAY THE TIMER
    if verbose:
        start_display = 'Process started at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")

    #ASSIGN THE FASTA FILE TO AN OBJECT
    aln = []

    with open(infile, 'r') as fastaFile:
        for line in fastaFile:
            #print line,
            aln.append(line)
    #print aln

    #REPLACE THE AMBIGUOUS CHARACTER STATES WITH N
    ambigs = 'rmwskyvhdb'.upper() # These are the ambiguous character states, from R package 'ape'
    ambigs = list(ambigs)
    missCount = [] # make an empty list to store the counts of missing data

    outfile = open(outfile, 'w') # open the outfile

    for line in tqdm(aln):
        if line.startswith('>'):
            #print line,
            outfile.write(line)
        else:
            seq = list(line)
            for nt in range(len(seq)):
                if seq[nt] in ambigs:
                    seq[nt] = 'N'
                else:
                    pass
            line = ''.join(seq)
            outfile.write(line)
            if verbose:
                #print line.count('N')
                C = line.count('N')
                missCount.append(C)
    #print numpy.percentile(missCount, 25)
    if verbose:
        print('\nMin N\t\t=\t{!s}'.format(numpy.min(missCount)))
        print('Q25 N\t\t=\t{!s}'.format(numpy.percentile(missCount, 25)))
        print('Mean N\t\t=\t{!s}'.format(numpy.mean(missCount)))
        print('Median N\t=\t{!s}'.format(numpy.median(missCount)))
        print('Q75 N\t\t=\t{!s}'.format(numpy.percentile(missCount, 75)))
        print('Max N\t\t=\t{!s}'.format(numpy.max(missCount)))
        time_adj = time.time() - start_time
        if time_adj < 60:
            print('*\tProcess completed in {!s} seconds'.format(
                round(time_adj, 2)))
        if time_adj >= 60 and time_adj < 3600:
            print('*\tProcess completed in {!s} minutes'.format(
                round((time_adj / 60), 2)))
        if time_adj >= 3600 and time_adj < 86400:
            print('*\tProcess completed in {!s} hours'.format(
                round(((time_adj / 60) / 60), 2)))
        if time_adj >= 86400:
            print('*\tProcess completed in {!s} days\n'.format(
                round((((time_adj / 60) / 60) / 24), 2)))
        print

if __name__=="__main__":
    main()

''' TEST PLESE IGNORE
aln = '>T1\nACGT\n>T2\nACGM\n>T3\nAYMT'
aln = aln.split()
print aln

#ambigs = ['Y', 'M']
ambigs = 'rmwskyvhdb'.upper()
ambigs = list(ambigs)

#outfile = '/Users/u5015730/Desktop/testRMambig.fasta'
#outfile = open(outfile, 'w')

for line in aln:
    if line.startswith('>'):
        print line
        #outfile.write(line + '\n')
    else:
        seq = list(line)
        for nt in range(len(seq)):
            if seq[nt] in ambigs:
                seq[nt] = 'N'
            else:
                pass
        line = ''.join(seq)
        print line
        #outfile.write(line + '\n')


outfile.close()


with open('/Users/u5015730/Desktop/testRMambig_orig.fasta', 'r') as F:
    for line in F:
        print line,
'''
