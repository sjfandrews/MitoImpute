#!/usr/bin/env python

import time
import os
import sys
import gzip
import argparse
import fcntl, termios, struct
import itertools
from tqdm import *

def memory_usage_resource():
    # FROM: http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
    import resource
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem

def terminal_size():
    # FROM USER pascal https://stackoverflow.com/users/362947/pascal
    # https://stackoverflow.com/questions/566746/how-to-get-linux-console-window-width-in-python
    h, w, hp, wp = struct.unpack('HHHH',
        fcntl.ioctl(0, termios.TIOCGWINSZ,
        struct.pack('HHHH', 0, 0, 0, 0)))
    return w, h

def main():
    
    start_time = time.time()
    start_display = 'Process started at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")
    
    #DEFINE INPUT ARGUMENTS
    parser=argparse.ArgumentParser(description='Curate gaps forced into the rCRS by newly aligned sequences',formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    
    parser.add_argument('-f', '--fasta', dest='fasta', type=str, required=True, help='specify the input FASTA file')
    parser.add_argument('-o', '--out', dest='out_dir', type=str, required=False, help='specify the output directory')
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", required=False, help='turn on verbose mode')
    
    args=parser.parse_args()
    
    fasta_file = args.fasta
    out_dir = args.out_dir
    verbose = args.verbose
    
    if verbose:
        print
        print(' START FASTA rCRS CURATION '.center(int(terminal_size()[0]), '='))
        print        
    
    if out_dir is None:
        out_dir = os.path.dirname(fasta_file)
        if verbose:
            print "*\tOUTPUT DIRECTORY NOT SPECIFIED"
            print "*\tDEFAULTING TO:\t%s" % (out_dir + "/")
            print
    if out_dir.endswith('/'):
        pass    
    else:
        out_dir += '/'
        
    outfile = os.path.basename(fasta_file)[:-6]
    outfile = out_dir + outfile + "_curated.fasta"
    
    #
    tkr = 0
    seq_names = []
    seqs = []
    rCRS_gaps = []
    rCRS_gaps2 = []
    n_correct = 0
    
    with open(fasta_file, 'r') as ff:
        for line in ff:
            tkr += 1
            if verbose:
                if tkr % 10000 == 0:
                    print "*\tPARSING FASTA LINE %s" % tkr
            if line.startswith(">"):
                line = line.strip("\n")
                seq_names.append(line)
            else:
                line = line.strip("\n")
                line = list(line)
                seqs.append(line)
            if tkr == 2:
                for n in range(len(line)):
                    if line[n] == "-":
                        rCRS_gaps.append(n)
    if verbose:
        print
        print "*\tCURATING OUT GAPS"
        for i in tqdm(range(len(seqs))):
            for j in rCRS_gaps:
                seqs[i][j] = ""
            seqs[i] = "".join(seqs[i])
            if len(seqs[i]) == 16569:
                n_correct += 1
            else:
                print "*\tERROR!"
                print "*\tERROR!"
                print "*\tERROR!"
                print "*\tSEQUENCE %s (%s) DID NOT CURATE TO CORRECT LENGTH" % (i ,seq_names[i])       
    else:
        for i in range(len(seqs)):
            for j in rCRS_gaps:
                seqs[i][j] = ""
            seqs[i] = "".join(seqs[i])
            if len(seqs[i]) == 16569:
                n_correct += 1
            else:
                print "*\tERROR!"
                print "*\tERROR!"
                print "*\tERROR!"
                print "*\tSEQUENCE %s (%s) DID NOT CURATE TO CORRECT LENGTH" % (i ,seq_names[i])                
            
    with open(outfile, 'wr') as of:
        if verbose:
            print "*\tWRITING TO FILE"
            for i in tqdm(range(len(seq_names))):
                of.write(seq_names[i] + "\n")
                of.write(seqs[i] + "\n")                
        else:
            for i in range(len(seq_names)):
                of.write(seq_names[i] + "\n")
                of.write(seqs[i] + "\n")    
    #
    
    stop_time = time.time()
    stop_display = 'Process completed at ' + time.strftime("%d/%m/%Y") + " " + time.strftime("%H:%M:%S")
    mem_used = memory_usage_resource()    
    
    if verbose:
        print
        print '*\tINPUT FASTA:\t\t\t%s' % fasta_file
        print '*\tOUTPUT FASTA:\t\t\t%s' % outfile
        print
        print "*\tNUMBER OF GAPS FORCED IN rCRS:\t%s" % len(rCRS_gaps)
        print '*\tNUMBER OF SAMPLES:\t\t%s' % len(seq_names)
        print "*\tNUMBER OF SAMPLES @ 16,569:\t%s" % n_correct
        print
        print '*\t' + str(start_display)
        print '*\t' + str(stop_display)
        time_adj = time.time() - start_time
        if time_adj < 60:
            print('*\tProcess completed in %s seconds' % (round(time_adj, 2)))
        if time_adj >= 60 and time_adj < 3600:
            print('*\tProcess completed in %s minutes' % (round((time_adj / 60), 2)))
        if time_adj >= 3600 and time_adj < 86400:
            print('*\tProcess completed in %s hours' % (round(((time_adj / 60) / 60), 2)))
        if time_adj >= 86400:
            print('*\tProcess completed in %s days' % (round((((time_adj / 60) / 60) / 24), 2)))
        if mem_used < 1024:
            print '*\tProcess used %s MB of memory' % ('%.2f' % (mem_used))
        if mem_used >= 1024:
            print '*\tProcess used %s GB of memory' % ('%.2f' % (mem_used / 1024))
        print
        print(' FINISH FASTA rCRS CURATION '.center(int(terminal_size()[0]), '='))
        print
    else:
        print
        print '*\tINPUT FASTA:\t\t%s' % fasta_file
        print '*\tOUTPUT FASTA:\t\t%s' % outfile
        print
        print '*\tNUMBER OF SAMPLES:\t%s' % len(seq_names)
        print
        print '*\t' + str(start_display)
        print '*\t' + str(stop_display)
        time_adj = time.time() - start_time
        if time_adj < 60:
            print('*\tProcess completed in %s seconds' % (round(time_adj, 2)))
        if time_adj >= 60 and time_adj < 3600:
            print('*\tProcess completed in %s minutes' % (round((time_adj / 60), 2)))
        if time_adj >= 3600 and time_adj < 86400:
            print('*\tProcess completed in %s hours' % (round(((time_adj / 60) / 60), 2)))
        if time_adj >= 86400:
            print('*\tProcess completed in %s days' % (round((((time_adj / 60) / 60) / 24), 2)))
        if mem_used < 1024:
            print '*\tProcess used %s MB of memory' % ('%.2f' % (mem_used))
        if mem_used >= 1024:
            print '*\tProcess used %s GB of memory' % ('%.2f' % (mem_used / 1024))
        print        

if __name__=="__main__":
    main()    

'''

infile = "/Volumes/TimMcInerney/SANDBOX/0_McInerney_Master_Alignment_July18_2018_1.fasta"
outfile = infile[:-6] + "curated.fasta"

tkr = 0
seq_names = []
seqs = []
rCRS_gaps = []
rCRS_gaps2 = []

with open(infile, 'r') as ff:
    for line in ff:
        tkr += 1
        if line.startswith(">"):
            line = line.strip("\n")
            seq_names.append(line)
        else:
            line = line.strip("\n")
            line = list(line)
            seqs.append(line)
        if tkr == 2:
            for n in range(len(line)):
                if line[n] == "-":
                    rCRS_gaps.append(n)

for i in range(len(seqs)):
    for j in rCRS_gaps:
        seqs[i][j] = ""
    seqs[i] = "".join(seqs[i])
    if len(seqs[i]) != 16569:
        print False
        
with open(outfile, 'wr') as of:
    for i in range(len(seq_names)):
        of.write(seq_names[i] + "\n")
        of.write(seqs[i] + "\n")
'''