#!/usr/bin/env python

import sys

guides_filename = sys.argv[1]
Fprimers_filename = sys.argv[2]
structuralRprimer_filename = sys.argv[3]
retron_filename = sys.argv[4]

with open(structuralRprimer_filename, 'r') as infile:
    structuralRprimer = infile.read().strip()

with open(retron_filename, 'r') as infile:
    retron = infile.read().strip()

# currently only works with 2 sublibraries
with open(Fprimers_filename, 'r') as infile:
    Fprimers = infile.read().split("\n")[1:4:2]

with open(guides_filename, 'r') as infile:
    for line in infile:
        if line.startswith("ID"):
            continue
        lineText = line.strip().split()
        oligoName = lineText[0]
        if oligoName.startswith("pan"):
            Fprimer = Fprimers[1]
        else:
            Fprimer = Fprimers[0]
        guide = lineText[10][:20]
        repairTemplate = lineText[13]
        full_oligo = Fprimer + repairTemplate + retron + guide + structuralRprimer
        print(oligoName + "\t" + full_oligo)