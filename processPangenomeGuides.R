#!/usr/bin/env Rscript
# R 3.6.3

library(data.table)

guideScores <- fread('data/pangenomeguides.scores.pam.tsv')
guideTargetCounts <- fread('data/guideTargetCounts.tsv')
setnames(guideTargetCounts, c("guide", "coreGenomePlus", "coreGenomeMinus", "panGenomePlus", "panGenomeMinus"))

threePrimeRetron <- "AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCA"
structuralgRNAprimer <- "GTTTCAGAGCTATGCTGGAA"



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

# truncate percent into cds to 2 decimal points
guides[, pctIntoCDS := floor(pctIntoCDS*100)/100]

# calculate "meta score" equal to (truncated sgRNAScorer2 value) * (1 - pctIntoCDS) * (10*PAM_count)

guides[, metaScore := 0]

# Prefer guides with score >= 0
guides[Score >= 0, metaScore := metaScore + 100000]

# Prefer guides with fewer PAMs in the sequence
guides[,metaScore := metaScore + ((10 - PAM_count) * 10000)]

# Prefer guides within first 50% of coding sequence
guides[pctIntoCDS <= 0.5, metaScore := metaScore + 1000]

# Break ties based on guide score bin
guides[,guideScoreBin := cut(Score, breaks=999)]
guides[,guideScoreBin := factor(guideScoreBin)]
guides[,guideScoreBin := as.numeric(guideScoreBin)]
guides[, metaScore := metaScore + guideScoreBin]

# rank each gene's gRNAs by metaScore
guides[, "metaRank" := frank(-metaScore), by=geneHeader]

# take all guides with rank <= 5
guides.final <- guides[metaRank <= 5]

# Check how many genes are covered
guides.final[, .N, by=geneHeader]

# Clean up table
guides.final[, c("SeqID", "PAM_count", "metaScore", "guideScoreBin", "metaRank") := NULL]

# Write output
fwrite(guides.final, "data/guides.final.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")