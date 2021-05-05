#!/usr/bin/env Rscript

# ensure R/3.6 and python/2.7 modules are loaded

library(data.table)
library(ggplot2)
library(ggthemes)


if(! file.exists("data/coreguides.tsv.gz")) {
    # load core guides
    coreGuides <- fread('zcat data/1011Genomes.guides.tsv.gz')

    # load gene lengths .tsv and reformat
    gene_lengths <- fread('data/gene_lengths.tsv')
    setnames(gene_lengths, c("ID2", "GENEID", "AA_length"))
    gene_lengths[, "ID2" := NULL]

    # merge gene AA lengths into guide table
    setkey(gene_lengths, GENEID)
    setkey(coreGuides, GENEID)
    coreGuides <- merge(coreGuides, gene_lengths)
    # Note that AA length is simply equal to codons.from.start + codons.from.stop

    # calculate how far into gene the guide targets
    coreGuides[, fraction.into.gene := codons.from.start / AA_length]

    # add unique ID for each guide
    coreGuides[, ID := 1:.N]

    # Score core genome guides (if scores do not yet exist)
    if(! file.exists("data/coreguides.scores.tab")) {


        # write guides (and their ID numbers) to .tsv file
        fwrite(coreGuides[, c("ID","guide")], file="coreguides.tsv", quote=F, row.names=F, col.names=F, sep="\t")

        # convert .tsv to .fasta
        system("sed 's/^/>/g' coreguides.tsv | tr '\t' '\n' > sgRNAScorer2/coreguides.fasta && rm coreguides.tsv")

        system("(cd sgRNAScorer2 && python identifyAndScore.py -i coreguides.fasta -o ../data/coreguides.scores.tab -p 3 -s 20 -l NGG)")   
        system("rm sgRNAScorer2/coreguides.fasta")
    }

    coreGuideScores <- fread('data/coreguides.scores.tab')

    # generate histogram of scores
    # hist(coreGuideScores$Score, breaks=50)

    # convert SeqID column back to simply numeric ID

    coreGuideScores[, ID := tstrsplit(SeqID, "_")[1]]
    coreGuideScores[, ID := as.numeric(ID)]
    coreGuideScores[, c("SeqID","Sequence") := NULL]
    setkey(coreGuideScores, ID)

    # Merge score into big table
    setkey(coreGuides, ID)
    coreGuides <- merge(coreGuideScores, coreGuides)

    # Write final guide score table to file
    fwrite(coreGuides, file="data/coreguides.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    system("gzip data/coreguides.tsv")
} else {
    coreGuides <- fread("zcat data/coreguides.tsv.gz")
}


# Flag putative guides where fraction.into.gene < 0.8 and sgRNA score > 0

dat[, QC_PASS := ifelse(fraction.into.gene < 0.8 & Score > 0 & frac.perfect > 0.99, TRUE, FALSE)]

# 477682 guides in first 80% of gene
# 163433 guides have score > 0
# 130419 guides pass both QC

# Check number of guides that pass both QC for various degrees of desired frac.perfect
o <- foreach(i = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.999, 1), .combine="rbind") %do% {
    data.table("frac.perfect"=i, "N"=nrow(dat[fraction.into.gene < 0.8 & Score > 0 & frac.perfect >= i]))
}

# Look at distribution for # of guides for each gene
ggplot(o, aes(x=frac.perfect, y=N)) +
geom_point() +
geom_label(aes(y=N+4000,x=frac.perfect+0.02,label=frac.perfect)) +
labs(x="fraction of strains covered by guide exactly", y="Number of guides")


dat[fraction.into.gene < 0.8 & Score > 0 & frac.perfect >= 0.8, .N, by=list(GENEID]

o <- foreach(i = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.999, 1), .combine="rbind") %do% {
    foreach(MinGuidesPerGene = 1:10, .combine="rbind") %do% {
        data.table("frac.perfect"=i, 
                    "MinGuidesPerGene" = MinGuidesPerGene,
                    "N_Genes"=nrow(dat[fraction.into.gene < 0.8 & Score > 0 & frac.perfect >= i, .N, by=list(GENEID)][N >= MinGuidesPerGene]))
    }
}


ggplot(o, aes(x=MinGuidesPerGene, y=N_Genes)) +
geom_point() +
facet_wrap(.~frac.perfect, labeller="label_both", nrow=2) +
labs(x="Minimum desired number of quality guides per gene", y="Number of genes matching all QC criteria") +
scale_x_continuous(breaks=1:10) +
theme_few(15)

# Try top 20 genes within frac.perfect first

dat.filtered <- dat[QC_PASS == TRUE]
dat.filtered[, "frac.perfect.rank" := frank(-frac.perfect, ties.method="random"), by=GENEID]
dat.filtered[frac.perfect.rank <= 20][, .N, by=GENEID][N>=5]

ggplot(dat.filtered[frac.perfect.rank <= 20], aes(x=GENEID, y=frac.perfect)) + geom_point()
