# lewis-collab

## download data
```
./getData.sh
```


## Build list of core gene IDs in our core gene guide table
```
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
```
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
```
python findGuidesFromFasta.py data/allORFS_pangenome.fasta data/pangenomeSpecificIDs.txt | \
    gzip > data/pangenome.guides.tsv.gz
```

## score pangenome guides
```
module purge
module load biopython
zcat data/pangenome.guides.tsv.gz | \
    awk 'NR > 1 {print $9,$5}' | \
sed 's/^/>/g' | tr ' ' '\n' > sgRNAScorer2/pangenomeguides.fasta

(cd sgRNAScorer2 && \
    python identifyAndScore.py \
    -i pangenomeguides.fasta \
    -o ../data/pangenomeguides.scores.tab \
    -p 3 \
    -s 20 \
    -l NGG)
```

## Process scored guides in R

```

```



## Structural reverse primer sequence
`GTTTCAGAGCTATGCTGGAA` sequence stored as file `data/structuralRprimer.dna`

## 3' Retron sequence
`AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCA` sequence stored as file `data/retron.dna`

