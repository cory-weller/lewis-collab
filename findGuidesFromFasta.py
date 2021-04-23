#!/usr/bin/env python

# Note: requires the `regex` module, which allows for overlapping matches
# bioconda/py3.7
# pip install regex

import regex as re
import sys

preNucleotides = 52
armLength = 50
guideLength = 20
pattern = re.compile("""[ACTG]{%s}[ACTG]GG[ACTG]{%s}""" % (preNucleotides, preNucleotides))
guideStartIdx = preNucleotides - guideLength
guideEndIdx = preNucleotides

inFileName = sys.argv[1]
desiredIDsFilename = sys.argv[2]


# load fasta

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
    }
    return(''.join(rev[x] for x in dnaSeq[::-1]))

def splitFasta(filename):
    n = 0
    with open(filename, "r") as infile:
        records = infile.read().split(">")[1:]
    for i in records:
        yield(fasta(i, n))
        n += 1

class fasta:
    help = 'stores type, header, and sequence information for FASTA files'
    def __init__(self, record, n):
        record = record.split("\n")
        self.header = record[0]
        # N is the number contained within the header, which is not perfectly sequential
        self.N = record[0].split("-")[0]
        # n is the index of the entry, starting at 0
        self.n = n
        self.seq = ''.join(x.rstrip() for x in record[1:])
        self.revCompSeq = revComp(self.seq)
        self.ntLength = len(self.seq)
        self.aaLength = self.ntLength / 3.0
        plusMatches = list(re.finditer(pattern, self.seq, overlapped=True))
        minusMatches = list(re.finditer(pattern, self.revCompSeq, overlapped=True))
        self.guides = [guide(x, "+") for x in plusMatches] + [guide(x, "-") for x in minusMatches]

class guide:
    help = 'stores index, strand, and repair template of guide RNA'
    def __init__(self, match, strand):
        self.PAMstrand = strand
        self.start = match.span()[0]
        self.end = match.span()[1]
        self.gRNA = match.group()[guideStartIdx:guideEndIdx]
        self.frame = (preNucleotides + self.start)%3
        if self.PAMstrand == "+" :
            self.repairTemplate = match.group()[((preNucleotides - armLength) - self.frame) : (preNucleotides - self.frame)].lower() + "TAA" + match.group()[(preNucleotides + 4 - self.frame):(preNucleotides + 4 - self.frame + armLength)].lower()
        elif self.PAMstrand == "-" :
            if self.frame == 0:
                self.offset = 1
            elif self.frame == 1:
                self.offset = 2
            elif self.frame == 2:
                self.offset = 0
            self.repairTemplate = match.group()[ (preNucleotides - armLength - self.offset) : (preNucleotides - self.offset) ].lower() + "TTA" + match.group()[(preNucleotides + 4 - self.offset) : (preNucleotides +4 - self.offset + armLength) ].lower()



fastaRecords = splitFasta(inFileName)
fastaRecords = list(fastaRecords)

# load desired IDs
with open(desiredIDsFilename, "r") as infile:
    desiredIDs = [x.strip() for x in infile.readlines()]

# print header
headerText = "\t".join([
    "geneHeader",
    "IDnumber",
    "ntLength",
    "aalength",
    "guide",
    "PAMstrand",
    "frame",
    "repairTemplate",
    "pctIntoCDS",
    "uniqueID"
])

print(headerText)
uniqueID = 1
# iterate through genes
for gene in fastaRecords:
    if (gene.N not in desiredIDs) or (gene.ntLength %3 != 0):
        continue
    # iterate through guides
    headerWithoutN = ''.join(gene.header.split("-")[1:])
    for guide in gene.guides:
        pctIntoGene = (guide.start + 50) / float(gene.ntLength)
        if (float(pctIntoGene) > 0.8):
            continue
        output = "\t".join([str(x) for x in [
            headerWithoutN,
            gene.N,
            gene.ntLength,
            int(gene.aaLength),
            guide.gRNA,
            guide.PAMstrand,
            guide.frame,
            guide.repairTemplate,
            pctIntoGene,
            uniqueID
            ]])
        print(output)
        uniqueID += 1


# fastaRecords[0].plusMatches[0]
# Print records with nt length indivisible by 3
# Likely to have introns, be pseudogenes, rRNA, transposons, etc
# a = list(fastaRecords)
# for gene in a:
#     if gene.ntLength%3 != 0:
#         print(gene.header.split("-")[1])

# Mutation Table
# 123|NGG|456 -> 123|TAA|-56
# 12N|GG3|456 -> 12N|TAA|-56
# 1NG|G23|456 -> 1NG|TAA|-56

# for gene in fastaRecords[0:4]:
#     for guide in gene.guidesF:
#         start = guide.span()[0]
#         end = guide.span()[1]
#         guideSeq = guide.group()[31:54]
#         offset = start%3
#         if offset == 0:
#             stopIndex = 51
#         elif offset == 1:
#             stopIndex = 50
#         elif offset == 2:
#             stopIndex = 52
#         repairTemplate = guide.group()[0:stopIndex].lower() + "TAA" + guide.group()[(stopIndex+4):].lower()
#         # If offset = 0, introduce STOP at 53, 54, 55
#         # If offset = 1, STOP at 5
#         # homology_right =
#         # Note: add 1 to index, as python is 0-indexed but we want 1-indexed
#         # output = "\t".join([str(x) for x in [gene.N, guide.span()[0]+1, guide.span()[1]+1, guide.group()]])
#         output = "\t".join([str(x) for x in [gene.N, guideSeq, repairTemplate]])
#         print(output)

# 180 / 7796 are not divisible by 3


