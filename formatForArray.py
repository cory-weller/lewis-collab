#!/usr/bin/env python

import sys

guides_filename = sys.argv[1]
Fprimer_filename = sys.argv[2]
structuralRprimer_filename = sys.argv[3]
retron_filename = sys.argv[4]

with open(structuralRprimer_filename, 'r') as infile:
    structuralRprimer = infile.read().strip()

with open(retron_filename, 'r') as infile:
    retron = infile.read().strip()

with open(Fprimer_filename, 'r') as infile:
    Fprimer = infile.read().strip()

with open(guides_filename, 'r') as infile:
    for line in infile:
        if line.startswith("uniqueID"):
            continue
        lineText = line.strip().split()
        oligoName = "lewis_" + lineText[0]
        guide = lineText[6][:20]
        repairTemplate = lineText[11]
        full_oligo = Fprimer + repairTemplate + retron + guide + structuralRprimer
        print(oligoName + "\t" + full_oligo)