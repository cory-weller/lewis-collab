#!/usr/bin/env python

# This script counts the number of PAM
import regex as re
import sys

data_infile = sys.argv[1]
# add 3 to account for PAM
nt_overlap_threshold = int(sys.argv[2]) + 3

PAM_plus = re.compile(r'[ACTG]GG')
PAM_minus = re.compile(r'CC[ACTG]')

with open(data_infile, 'r') as infile:
    for line in infile:
        if line.startswith("SeqID"):
            print("SeqID\tSequence\tScore\tPAM_count")
            continue
        SeqID, Sequence, Score = line.split()
        N_plus = len(re.findall(PAM_plus, Sequence[-nt_overlap_threshold:], overlapped = True))
        N_minus = len(re.findall(PAM_minus, Sequence[:nt_overlap_threshold], overlapped = True))
        N = N_plus + N_minus
        print('\t'.join([str(x) for x in [SeqID, Sequence, Score, N]]))


