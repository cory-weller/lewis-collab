#!/usr/bin/env python3

import sys
import random

random.seed(1)

gRNAlength = 20
RTlength = 103
n = 100

retron = "AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCA"
structuralRprimer = "GTTTCAGAGCTATGCTGGAA"
fPrimer = "CGATCGCCCTTGGTG"

nucleotides = ['A','C','T','G']

def wrapFasta(seq):
    return '\n'.join([seq[x:x+80] for x in range(0,len(seq),80)])


for i in range(1, n+1):
    RT = ''.join(random.choices(nucleotides, k=RTlength))
    gRNA = ''.join(random.choices(nucleotides, k=gRNAlength))
    oligo = fPrimer + RT + retron + gRNA + structuralRprimer
    print("control_oligo-%s\t%s" % (i, oligo))