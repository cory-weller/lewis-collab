# README

## download data
```bash
./getData.sh
```

## Process core guides
```bash
module load R/3.6.3
Rscript processCoreGuides.R
```

## Build list of core gene IDs in our core gene guide table
```bash
zcat data/1011Genomes.guides.tsv.gz | \
    awk 'NR > 1 {print $8}' | \
    sort -u > data/1011Genomes.coreGenes.txt

# check number of genes in core genome set:
# wc -l data/1011Genomes.coreGenes.txt
# 6488
```

## Find sequences unique to pangenome
Inverse match core gene IDs on pangenome IDs
Total number of genes generated based on core genome:
```bash
grep ">" data/allORFS_pangenome.fasta | \
    grep -v -f data/1011Genomes.coreGenes.txt | \
    cut -d ">" -f 2 > data/pangenomeSpecificGenes.txt

# check number of genes unique to this set:
# wc -l data/pangenomeSpecificGenes.txt
# 1783

# save list of numeric IDs for this pangenome-specific genes
cut -d "-" -f 1 data/pangenomeSpecificGenes.txt > data/pangenomeSpecificIDs.txt
```

## Generate pangenome-specific guide sequences
```bash
module purge python
module load python/3.6
python findGuidesFromFasta.py data/allORFS_pangenome.fasta data/pangenomeSpecificIDs.txt | \
    gzip > data/pangenome.guides.tsv.gz
```

## Score pangenome guides
```bash
module purge python
module load python/2.7
# awk column 12 is unique ID, column 5 is guide sequence
zcat data/pangenome.guides.tsv.gz | \
    awk 'NR > 1 {print $12,$5}' | \
sed 's/^/>/g' | tr ' ' '\n' > sgRNAScorer2/pangenomeguides.fasta


# Note sgRNAScorer2 requires Python 2.7 or greater; python packages: BioPython and and scikit-learn
# Note guide RNA sequences must include the NGG PAM sequence
module purge python
module load python/2.7
(cd sgRNAScorer2 && \
    python identifyAndScore.py \
    -i pangenomeguides.fasta \
    -o ../data/pangenomeguides.scores.tab \
    -p 3 \
    -s 20 \
    -l NGG)
```

## Retrieve reference sequences from NCBI
```bash
# Download S288c reference genome fasta files:
module load edirect
efetch -db nucleotide -id KP263414 -format fasta  > data/S288c.fasta
efetch -db nucleotide -id BK006935 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006936 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006937 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006938 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006939 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006940 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006941 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006934 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006942 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006943 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006944 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006945 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006946 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006947 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006948 -format fasta >> data/S288c.fasta 
efetch -db nucleotide -id BK006949 -format fasta >> data/S288c.fasta 

# BK006935 I
# BK006936 II
# BK006937 III
# BK006938 IV
# BK006939 V
# BK006940 VI
# BK006941 VII
# BK006934 VIII
# BK006942 IX
# BK006943 X
# BK006944 XI
# BK006945 XII
# BK006946 XIII
# BK006947 XIV
# BK006948 XV
# BK006949 XVI
# KP263414 mitochondrion
```

## Calculate multiple-targetting of pangenome guides
```bash
sbatch --time=12:00:00	guideTargeting.sh
```

## Count number of PAM matches in pangenome guides
This is to avoid using guides that have too much overlapping sequence (10 bp or more)
```bash
python3 countPAM.py "data/pangenomeguides.scores.tab" 10 > data/pangenomeguides.scores.pam.tsv
```

## Structural reverse primer sequence
`GTTTCAGAGCTATGCTGGAA` sequence stored as file `data/structuralRprimer.dna`

## 3' Retron sequence
`AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCA` sequence stored as file `data/retron.dna`

## Format oligos for array ordering
Note: Need to possibly change 5' primer sequence.

## Process pangenome guides and combine with core guides

```bash
module load R/3.6.3
Rscript processPangenomeGuides.R
```



```bash
python3 formatForArray.py \
    data/guides.merged.filter.tsv \
    data/Fprimers.dna \
    data/structuralRprimer.dna \
    data/retron.dna \
    > yeast-ko-array.txt
```

## Add control oligos
```bash
python3 addControlOligos.py >> yeast-ko-array.txt
```
