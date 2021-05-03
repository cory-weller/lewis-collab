#!/usr/bin/env Rscript
# R 3.6.3

library(data.table)

guideScores <- fread('data/pangenomeguides.scores.tab')
guideTargetCounts <- fread('data/guideTargetCounts.tsv')
setnames(guideTargetCounts, c("guide", "coreGenomePlus", "coreGenomeMinus", "panGenomePlus", "panGenomeMinus"))

panGenomeGuides <- cbind(guideTargetCounts, guideScores)
panGenomeGuides[, guide := NULL]
panGenomeGuides[, Sequence := NULL]
panGenomeGuides[, uniqueID := tstrsplit(SeqID, "_")[1]]
panGenomeGuides[, uniqueID := as.numeric(uniqueID)]

RTs <- fread('zcat data/pangenome.guides.tsv.gz')

setkey(RTs, uniqueID)
setkey(panGenomeGuides, uniqueID)

guides.merged <- merge(panGenomeGuides, RTs)

# subset guides that do not hit core genome, and only hit pangenome once
guides.merged[, panGenomeTotal := panGenomePlus + panGenomeMinus]
guides <- guides.merged[coreGenomeMinus == 0 & coreGenomePlus == 0 & panGenomeTotal == 1]
guides[, c("panGenomeTotal","coreGenomePlus","coreGenomeMinus","panGenomePlus","panGenomeMinus") := NULL]
# left with 74138 guides for pangenome

setkey(guides, endGuide)

# exclude guides ending in GGG, which 
guides <- guides[ ! grepl("GGG$", guide, perl=TRUE)]

# calculate "meta score" equal to (sgRNAScorer2 value) * (1 - pctIntoCDS)
guides[, metaScore := Score * (1-pctIntoCDS)]

# rank each gene's gRNAs by metaScore
guides[, "metaRank" := frank(-metaScore), by=geneHeader]

# take all guides with rank <= 5
guides.final <- guides[metaRank <= 5]

# Check how many genes are covered
guides.final[, .N, by=geneHeader]
# 996 genes have at least one guide


## HOW TO ADDRESS OVERLAP QUESTION?


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