#!/usr/bin/env Rscript
# R 3.6.3

library(data.table)

guideScores <- fread('data/pangenomeguides.scores.tab')
guideScores[, uniqueID := tstrsplit(SeqID, "_")[1]]
guideScores[, uniqueID := as.numeric(uniqueID)]
setkey(guideScores, uniqueID)

guides <- fread('zcat data/pangenome.guides.tsv.gz')

setkey(guides, uniqueID)

dat <- merge(guides, guideScores)

dat <- dat[Score > 0]

dat[, .N, by=list(geneHeader)]

atLeastFive <- dat[, .N, by=list(geneHeader)][N >= 5]$geneHeader

dat[, scoreRank := frank(-Score), by=geneHeader]

dat <- dat[geneHeader %in% atLeastFive & scoreRank <= 5]

dat [, c("SeqID","Sequence","scoreRank", "uniqueID") := NULL]

fwrite(dat, file="data/pangenomeguides.final.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

817 of 1077 pangenome-specific specific sequences have 5 guides within  80% of gene