#!/usr/bin/env Rscript
# R 3.6.3

library(data.table)

if(! file.exists("data/guides.merged.tsv.gz")) {
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

    panGenomeGuides <- merge(panGenomeGuides, RTs)
    panGenomeGuides[, "Genome" := "pan"]

    # subset guides that do not hit core genome, and only hit pangenome once
    panGenomeGuides[, panGenomeTotal := panGenomePlus + panGenomeMinus]
    panGenomeGuides <- panGenomeGuides[coreGenomeMinus == 0 & coreGenomePlus == 0 & panGenomeTotal == 1]
    panGenomeGuides[, c("panGenomeTotal","coreGenomePlus","coreGenomeMinus","panGenomePlus","panGenomeMinus") := NULL]

    # Change start and end to always reflect + strand
    panGenomeGuides[PAMstrand == "-", "tmp" := startGuide]
    panGenomeGuides[PAMstrand == "-", endGuide := startGuide]
    panGenomeGuides[PAMstrand == "-", startGuide := tmp]
    panGenomeGuides[, "tmp" := NULL]

    # update endGuide to include the PAM as well
    panGenomeGuides[, endGuide := endGuide + 3]

    # format pangenome guides for merging
    panGenomeGuides[, "ID" := paste0(Genome, uniqueID)]
    panGenomeGuides[, c("uniqueID", "IDnumber") := NULL]
    panGenomeGuides[, "chrom" := "pangenome"]
    setnames(panGenomeGuides, "guide", "guideSeq")
    setnames(panGenomeGuides, "geneHeader", "geneID")
    setnames(panGenomeGuides, "aalength", "aaLength")
    panGenomeGuides[, "ntLength" := NULL]
    panGenomeGuides[, "frame" := NULL]
    setnames(panGenomeGuides, "Score", "guideScore")
    panGenomeGuides[, "frac.perfect" := "NA"]
    panGenomeGuides[, "min.intron.dist" := "NA"]
    panGenomeGuides[, "terminator" := "NA"]
    panGenomeGuides[, "SeqID" := NULL]

    mergeColOrder <-  c(
    "ID",
    "Genome",
    "chrom",
    "geneID",
    "aaLength",
    "terminator",
    "startGuide",
    "endGuide",
    "pctIntoCDS",
    "PAMstrand",
    "guideSeq",
    "guideScore",
    "PAM_count",
    "frac.perfect",
    "repairTemplate",
    "min.intron.dist"
    )


    setcolorder(panGenomeGuides, mergeColOrder)


    # Load core genome guides and format for merging
    coreGenomeGuides <- fread('zcat data/coreGuides_50bp.tsv.gz')
    coreGenomeGuides[, "Genome" := "core"]
    coreGenomeGuides[, "ID" := paste0(Genome, ID)]
    setnames(coreGenomeGuides, "guide", "guideSeq")
    setnames(coreGenomeGuides, "GENEID", "geneID")
    setnames(coreGenomeGuides, "AA_length", "aaLength")
    setnames(coreGenomeGuides, "Score", "guideScore")
    setnames(coreGenomeGuides, "fraction.into.gene", "pctIntoCDS")
    setnames(coreGenomeGuides, "repair.templates", "repairTemplate")
    setnames(coreGenomeGuides, "start", "startGuide")
    setnames(coreGenomeGuides, "end",   "endGuide")

    coreGenomeGuides[, c("CDSLOC.start", "CDSLOC.end", "codons.from.start", "codons.from.stop", "REFCODON", "REFAA" ) := NULL]
    coreGenomeGuides[, c("PROTEIN.start", "PROTEIN.end") := NULL]

    setcolorder(coreGenomeGuides, mergeColOrder)

    guides <- rbindlist(list(panGenomeGuides, coreGenomeGuides))

    fwrite(guides, file="data/guides.merged.prefilter.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    system("gzip data/guides.merged.prefilter.tsv")
} else {
    guides <- fread("zcat data/guides.merged.prefilter.tsv.gz")
}



# truncate percent into cds to 2 decimal points
guides[, pctIntoCDS := floor(pctIntoCDS*100)/100]

# calculate "meta score" equal to (truncated sgRNAScorer2 value) * (1 - pctIntoCDS) * (10*PAM_count)

guides[, metaScore := 0]

# Prefer guides with score >= 0
guides[guideScore >= 0, metaScore := metaScore + 100000]

# Prefer guides with fewer PAMs in the sequence
guides[,metaScore := metaScore + ((10 - PAM_count) * 10000)]

# Prefer guides within first 50% of coding sequence
guides[pctIntoCDS <= 0.5, metaScore := metaScore + 1000]

# Break ties based on guide score bin
guides[,guideScoreBin := cut(guideScore, breaks=999)]
guides[,guideScoreBin := factor(guideScoreBin)]
guides[,guideScoreBin := as.numeric(guideScoreBin)]
guides[, metaScore := metaScore + guideScoreBin]

# rank each gene's gRNAs by metaScore
guides[, "metaRank" := frank(-metaScore), by=list(Genome, geneID)]

# take all guides with rank <= 5
guides.final <- guides[metaRank <= 5]
guides.final[, ID := paste(Genome, geneID, metaRank, sep="-")]

# Check how many genes are covered
guides.final[, .N, by=geneID]

# Clean up table
guides.final[, c("PAM_count", "metaScore", "guideScoreBin", "metaRank") := NULL]

# Write output
fwrite(guides.final, "data/guides.merged.filter.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
