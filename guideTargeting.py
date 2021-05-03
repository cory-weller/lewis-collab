#!/usr/bin/env python

import sys
import regex as re

coreGenomeFilename = sys.argv[1]
panGenomeFilename = sys.argv[2]
guidesFilename = sys.argv[3]

# number of nucleotides (from 3' end) requiring perfect match
threePrimeNucleotides = 12

coreGenomeFilename = "data/S288c.fasta"
panGenomeFilename = "data/allORFS_pangenome.fasta"
guidesFilename = "data/pangenomeguides.scores.tab"

def concatenateSeqs(filename):
    with open(filename, 'r') as infile:
        # For pattern searching, concatenate all chromosomes, 
        # separated by 50 '-' characters
        seqs = infile.read().split('>')[1:]
        seqs = [''.join(x.strip().split("\n")[1:]) for x in seqs]
        seqs = (50 * '-').join(seqs)
        return(seqs)

def readGuides(filename):
    guides = []
    with open(filename, 'r') as infile:
        guides = [x.rstrip().split()[1][:-3] for x in  infile.readlines()[1:]]
    return(guides)


def revComp(dnaSeq):
    dnaSeq = dnaSeq.upper()
    rev = {
            "A" : "T",
            "T" : "A",
            "C" : "G",
            "G" : "C",
            "N" : "N",
            "Y" : "R",
            "R" : "Y",
            "W" : "W",
            "S" : "S",
            "K" : "M",
            "M" : "K",
            "D" : "H",
            "V" : "B",
            "H" : "D",
            "B" : "V",
            "X" : "X",
            "-" : "-",
    }
    return(''.join(rev[x] for x in dnaSeq[::-1]))

def countGuideMatches(guideSeq, genomeSeq):
    guidePattern = re.compile("%s[ACTG]GG" % guideSeq)
    matches = re.findall(guidePattern, genomeSeq, overlapped=True)
    return(len(matches))

coreGenome = concatenateSeqs(coreGenomeFilename)
coreGenomeRevComp = revComp(coreGenome)
panGenome = concatenateSeqs(panGenomeFilename)
panGenomeRevComp = revComp(panGenome)

guides = readGuides(guidesFilename)


for guide in guides:
    guideSeq = guide[(-1 * threePrimeNucleotides):]
    coreGenomeMatches = countGuideMatches(guideSeq, coreGenome)
    coreGenomeRevCompMatches = countGuideMatches(guideSeq, coreGenomeRevComp)
    panGenomeMatches = countGuideMatches(guideSeq, panGenome)
    panGenomeRevCompMatches = countGuideMatches(guideSeq, panGenomeRevComp)
    print("\t".join([str(x) for x in [
                                    guide, 
                                    coreGenomeMatches,
                                    coreGenomeRevCompMatches,
                                    panGenomeMatches,
                                    panGenomeRevCompMatches
                                    ]]))
