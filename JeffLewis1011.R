#strategy:

#1) Figure out which gRNAs will work in which strains

# Figure out the position of each gRNA
# Extract the part of the VCF corresponding to that gRNA
# For each individual, tell whether it is completely WT or not in that range

#2) Choose your gRNAs.

# Prioritize gRNAs that will work in many different strains, and that:
# lack a TTTT terminator
# Occur in the "middle" of genes (not too close to the start, and far from the end)
# Aren't close to other chosen gRNAs, in case a given genomic region has a problem.
# Aren't close to intron/exon boundaries (the repair template could be wrong for those where the PAM is split across an intron)

#3) Figure out the repair template you want to use for each gRNA

#Josh's table has each gRNA in its codon context. Make your mutation in this context.
#Then add on appropriately-sized arms to get sufficient homology.

#Mutation table:

# 123|NGG|456 -> 123|TAA|-56
# 12N|GG3|456 -> 12N|TAA|-56
# 1NG|G23|456 -> 1NG|TAA|-56
# CCN|123 -> TAA|-23
# 1CC|N23 -> TAA|-23
# 12C|CN3 -> TAA|-N3

#I downloaded Joseph's VCF file from http://1002genomes.u-strasbg.fr/files/ (1011Matrix.gvcf.gz)

#I locally installed homebrew, following https://github.com/Homebrew/brew/blob/master/docs/Installation.md#untar-anywhere:
#(The local install is only because I don't have administrator privileges on this computer; the normal homebrew installation would be much better)

# mkdir homebrew && curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C homebrew

#I then used the local homebrew to locally install samtools and vcftools

# cd ~/Downloads/homebrew/bin
# ./brew install samtools
# ./brew install vcftools

#I created a tabix index of Joseph's vcf

# ~/Downloads/homebrew/Cellar/htslib/1.9/bin/tabix -p ~/Downloads/1011Matrix.gvcf.gz

#Then, I load Joshs' database of gRNAs.

library(GenomicRanges)
library(S4Vectors) #For help installing these libraries, see the email titled "administrator privileges" on 10/8/18.
library(Biostrings)
load("~/data/others_data/josh/BYGuideTable.RData") #Downloaded from https://www.dropbox.com/s/vfyqgmtp77iz4fo/BYGuideTable.RData?dl=0

#Some processing to make a table of gRNAs that target genes
codinggRNAs <- rbind(as.data.frame(BYguideList$codingEffectsNoCs)[,-10], as.data.frame(BYguideList$codingEffectsCs)[,-10])
codinggRNAs <- codinggRNAs[,c(1:3, 6, 7, 10, 12, 13, 15, 19, 21, 23)]
codinggRNAs <- codinggRNAs[order(codinggRNAs$seqnames, codinggRNAs$start),]
codinggRNAs <- cbind.data.frame(chrom = rep(1:16, rle(as.vector(codinggRNAs$seqnames))$lengths), codinggRNAs[,-1]) #I checked first that the chromosomes do in fact appear in the correct order
uniquecodinggRNAs <- codinggRNAs[which(codinggRNAs$guidePAM12PMcount == 1),]
uniquecodinggRNAs[,-6] -> uniquegRNAs



###########
#This is just for reference of how I built the table of gRNA features - if you don't need to redo 
#it, skip to 'load("~/code/JeffLewis1011/uniquegRNAs.Rdata")' (currently line 276)
###########

#The function below will generate a mini-vcf from Joseph's master vcf that focuses on the 
#guide targeting and PAM region for a specific gRNA, and then see which strains are WT in that
#mini-vcf. You pass the appropriate features as so:

#You can get chrom, grna.start, grna.stop, and strand from the BYguideList$codingEffectsCs and
#BYguideList$codingEffectsNoCs structures, from seqnames, the start and end of ranges, and
#PAMstrand, respectively.

wt.check <- function (chrom, grna.start, grna.stop, strand) {
  if (strand == "+") { 
    grna.stop2 <- grna.start + 19
    grna.start2 <- grna.start + 21
  } else {
    grna.stop2 <- grna.start + 1
    grna.start2 <- grna.start + 3
  }
  system(paste(sep = "", 
               collapse = "", 
               "~/Downloads/homebrew/Cellar/htslib/1.9/bin/tabix ~/Downloads/1011Matrix.gvcf.gz chromosome",
               chrom,
               ":",
               grna.start,
               "-",
               grna.stop2,
               " chromosome",
               chrom,
               ":",
               grna.start2,
               "-",
               grna.stop,
               " > ~/data/others_data/joseph_1011/test4.vcf"))
  if (file.size("~/data/others_data/joseph_1011/test4.vcf")) {
    dump <- scan("~/data/others_data/joseph_1011/test4.vcf", what = "character")
    if (length(dump) > 1020) {
      grna.vcf <- as.data.frame(matrix(dump, nrow = length(dump)/1020, byrow = T)[,-c(1:9)], stringsAsFactors = F) #I decided these two lines were faster than read.table.
      apply(grna.vcf, 2, function (x) min(sapply(x, function (y) substr(y, 1, 3) == "0/0")))
    } else {
      as.numeric(substr(dump[-c(1:9)], 1, 3) == "0/0")
    }
  } else rep(1, 1011)
}

#Make the giant table

grna.perfect.in.strain <- matrix(nrow = 605799, ncol = 1011)

for (x in 1:605799) {
  grna.perfect.in.strain[x,] <- wt.check(uniquegRNAs$chrom[x],
                                         uniquegRNAs$start[x],
                                         uniquegRNAs$end[x],
                                         uniquegRNAs$PAMstrand[x])
}

colnames(grna.perfect.in.strain) <- strsplit(system("head -n 59 ~/data/others_data/joseph_1011/1011Matrix.gvcf | tail -n 1", intern = T), "\t")[[1]][10:1020]

frac.perfect <- rowSums(grna.perfect.in.strain)/ncol(grna.perfect.in.strain)

uniquegRNAs <- cbind.data.frame(uniquegRNAs, frac.perfect)

cbind.data.frame(uniquegRNAs, frac.perfect, grna.perfect.in.strain) -> grna.perfect.in.strain.table
save(grna.perfect.in.strain.table, file = "~/code/JeffLewis1011/grna.perfect.in.strain.table.2.Rdata")

load("~/code/JeffLewis1011/grna.perfect.in.strain.table.2.Rdata")

frac.perfect <- rowSums(grna.perfect.in.strain.table[,12:1022])/ncol(grna.perfect.in.strain.table[,12:1022])

#I made the protein length table below from yeast mine's query builder,
#selecting Protein, then on the next page selecting (under the "Genes"
#expandable selection menu) Standard Name and Systematic name, and under the
#"Protein" menu selecting Length (in that order), then exporting as tsv.

read.table("~/data/others_data/sgd/gene_lengths.tsv") -> gene_lengths
colnames(gene_lengths) <- c("GENE", "ORF", "Protein.length")

#Figuring out where the PAMs fall in the ORF (I use the 2nd G of the PAM as its coordinate,
#because that one is always in the changed codon ccording to the mutation table at the top of this file)

# If PAMstrand is '+' and gene is W, or '-' and C, then PAM is CDSLOC.end, otherwise it's CDSLOC.start 

guide.gene.orientation <- numeric(nrow(uniquegRNAs))
for (x in 1:nrow(uniquegRNAs)) if ((substr(uniquegRNAs$GENEID[x], 7, 7) == "W") == (uniquegRNAs$PAMstrand[x] == "+")) guide.gene.orientation[x] <- uniquegRNAs$CDSLOC.end[x] else guide.gene.orientation[x] <- uniquegRNAs$CDSLOC.start[x]

#Then based on gene_lengths, figure out how far the cut and PAM are from the 3' end of the gene.

uniquegRNAs <- cbind.data.frame(uniquegRNAs, codons.from.start = ceiling(guide.gene.orientation/3))

uniquegRNAs <- cbind.data.frame(uniquegRNAs, codons.from.stop = apply(uniquegRNAs, 1, function (x) gene_lengths$Protein.length[which(gene_lengths$ORF == x$GENEID)] - x$codons.from.start))

#Figuring out the distance from each PAM to the nearest intron/exon boundary

#I made the intron position table below from yeast mine's query builder,
#selecting Intron, then on the next page selecting (under the "Genes"
#expandable selection menu) Standard Name and Systematic name, and under the
#"Locations" menu selecting End and Start, then exporting as tsv.

read.table("~/data/others_data/sgd/intron_positions.tsv") -> intron_positions
colnames(intron_positions) <- c("ORF", "GENE", "intron.end", "intron.start")

min.intron.dist <- apply(uniquegRNAs, 1, function (x) if (x$GENEID %in% intron_positions$ORF) {
  unlist(intron_positions[which(intron_positions$ORF == x$GENEID),c(3,4)]) -> intron.bounds
  if (x$PAMstrand == "+") pos <- x$end else pos <- x$start
  min(abs(pos - intron.bounds))
  } else Inf
  )

uniquegRNAs <- cbind.data.frame(uniquegRNAs, min.intron.dist)

#The repair template for each gRNA. I used the mutation table from lines 23-28.

#Because of problems installing bioconductor packages, I did this using seqinR.

#install.packages("seqinr")
library(seqinr)

#I downloaded the sacCer3 genome from https://hgdownload-test.gi.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
saccer3 <- list()
for (i in 1:16) saccer3[i] <- read.fasta(paste(sep = "", collapse = "", "~/data/others_data/sgd/sacCer3/", list.files("~/data/others_data/sgd/sacCer3/")[c(1:4, 7:10, 5, 11:17)[i]]))

#For gRNAs in 'W' genes, and with the gRNA on the '+' strand:

#If CDSLOC.end mod 3 is 0 (ie the NGG PAM is in frame), then the repair template is:
#end - 51 to end - 3, then TAA, then end + 2 to end + 50
#If CDSLOC.end %% 3 == 1, then the repair template is:
#end - 49 to end - 1, then TAA, then end + 4 to end + 52
#If CDSLOC.end %% 3 == 2, then the repair template is:
#end - 50 to end - 2, then TAA, then end + 3 to end + 51

#Shortcut to get there (for gRNAs in 'W' genes and with the gRNA on the '+' strand): 
#offset = ((CDSLOC.end + 2) %% 3) + 1
#end - (48 + offset) to end - offset, then TAA, then end + (5 - offset) to end + (53 - offset)

#For gRNAs in 'W' genes, and with the gRNA on the '-' strand:
#offset = ((CDSLOC.start + 2) %% 3) + 1
#start - (48 + offset) to start - offset, then TAA, then start + (5 - offset) to start + (53 - offset)

#For gRNAs in 'C' genes, and with the gRNA on the '-' strand:
#offset = ((CDSLOC.end + 2) %% 3) + 1
#start - (53 - offset) to start - (5 - offset), then TTA, then start + offset to start + 48 + offset

#For gRNAs in 'C' genes, and with the gRNA on the '+' strand:
#offset = ((CDSLOC.start + 2) %% 3) + 1
#end - (53 - offset) to end - (5 - offset), then TTA, then end + offset to end + 48 + offset

#So you get a function that looks like:
#(Test a bunch of repair templates once you've made them to make sure they're correct!!!!!)

reptemp <- function (gene.strand, PAM.strand, CDS.coord.start, CDS.coord.end, chrom.coord.start, chrom.coord.end, chrom) {
  if (gene.strand == "W" && PAM.strand == "+") {
    offset <- ((CDS.coord.end + 2) %% 3) + 1
    paste(paste(collapse = "",
                sep = "",
                saccer3[[chrom]][(chrom.coord.end - (48 + offset)):(chrom.coord.end - offset)]),
          "TAA",
          paste(collapse = "", 
                sep ="", 
                saccer3[[chrom]][(chrom.coord.end + (5 - offset)):(chrom.coord.end + (53 - offset))]), 
          collapse = "", sep ="")
  } else if (gene.strand == "W" && PAM.strand == "-") {
    offset <- ((CDS.coord.start + 2) %% 3) + 1
    paste(paste(collapse = "",
                sep = "",
                saccer3[[chrom]][(chrom.coord.start - (48 + offset)):(chrom.coord.start - offset)]),
          "TAA",
          paste(collapse = "", 
                sep ="", 
                saccer3[[chrom]][(chrom.coord.start + (5 - offset)):(chrom.coord.start + (53 - offset))]), 
          collapse = "", sep ="")
  } else if (gene.strand == "C" && PAM.strand == "+") {
    offset <- ((CDS.coord.start + 2) %% 3) + 1
    paste(paste(collapse = "",
                sep = "",
                saccer3[[chrom]][(chrom.coord.end - (53 - offset)):(chrom.coord.end - (5 - offset))]),
          "TTA",
          paste(collapse = "", 
                sep ="", 
                saccer3[[chrom]][(chrom.coord.end + offset):(chrom.coord.end + 48 + offset)]), 
          collapse = "", sep ="")
  } else if (gene.strand == "C" && PAM.strand == "-") {
    offset <- ((CDS.coord.end + 2) %% 3) + 1
    paste(paste(collapse = "",
                sep = "",
                saccer3[[chrom]][(chrom.coord.start - (53 - offset)):(chrom.coord.start - (5 - offset))]),
          "TTA",
          paste(collapse = "", 
                sep ="", 
                saccer3[[chrom]][(chrom.coord.start + offset):(chrom.coord.start + 48 + offset)]), 
          collapse = "", sep ="")
  }
}

repair.templates <- character(nrow(uniquegRNAs))

for (x in 1:nrow(uniquegRNAs)) repair.templates[x] <- reptemp(substr(uniquegRNAs$GENEID[x], 7, 7), 
                                                              uniquegRNAs$PAMstrand[x], 
                                                              uniquegRNAs$CDSLOC.start[x],  
                                                              uniquegRNAs$CDSLOC.end[x], 
                                                              uniquegRNAs$start[x], 
                                                              uniquegRNAs$end[x], 
                                                              uniquegRNAs$chrom[x])

uniquegRNAs <- cbind.data.frame(uniquegRNAs, frac.perfect)

uniquegRNAs <- cbind.data.frame(uniquegRNAs, terminator = grepl("TTTT", uniquegRNAs$guide))

uniquegRNAs <- cbind.data.frame(uniquegRNAs, repair.templates)

#save(uniquegRNAs, file = "~/code/JeffLewis1011/uniquegRNAs.Rdata")
#load("~/code/JeffLewis1011/uniquegRNAs.Rdata")

split(uniquegRNAs, uniquegRNAs$GENEID) -> uniquegRNAs.bygene

#Process for choosing gRNAs

#(1) Make a select list of gRNAs that will target >90% of strains, that do not have a terminator sequence, 
#that are more than 30 codons in from the start and 100 codons in from the end, and more than 5 bp away from a 
#splice junction.

#(2) For each gene, find the gRNA on that list that targets the largest number of strains; as a tiebreak,
#go with the gRNA with the smallest codons.from.start. Choose it.

#(3) Block out gRNAs that overlap in their targeting sequence with that gRNA. Reiterate. Go 5 times, or until
#you run out of gRNAs.

#(4) For genes without a complete set of 5 gRNAs, fill them in with similar criteria applied to a larger set of
#gRNAs: those >30 amino acids away from the end of the gene. Then add in gRNAs that overlap those that have
#already been chosen. After that, choose gRNAs outside the final 30 amino acids that are conserved in less than
#90% of strains, in order of decreasing conservation. After all that, if you still aren't at 5 gRNAs, add in gRNAs
#that are less than 30 amino acids from the end of the gene, in order of distance from the end of the gene.
#Note that there are 156 genes for which you don't have 5 gRNAs before any selection (other than the initial
#filter to remove gRNAs that target more than once).

selectgRNAs <- uniquegRNAs[(uniquegRNAs$terminator == F) & 
                             (uniquegRNAs$frac.perfect > .9) & 
                             (uniquegRNAs$codons.from.start > 29) & 
                             (uniquegRNAs$codons.from.stop > 99) & 
                             (uniquegRNAs$min.intron.dist > 4),]

split(selectgRNAs, selectgRNAs$GENEID) -> selectgRNAs.bygene

library(intervals)

pickgrnas <- function (gRNA.table) {
  gRNA.table <- gRNA.table[order(gRNA.table$codons.from.start),]
  preference <- order(gRNA.table$frac.perfect, decreasing = T)
  picked <- preference[1]
  rounds <- 1
  while (rounds < 5) { 
    preference <- preference[-which(preference %in% unique(unlist(interval_overlap(Intervals(gRNA.table[picked, c(2,3)]), Intervals(gRNA.table[,c(2,3)])))))] #this removes entries from consideration whose gRNAs overlap gRNAs already chosen
    if (length(preference)) { #Only proceed if there is anything left to consider
      picked <- c(picked, preference[1])
      rounds <- rounds + 1 
    } else rounds <- 5
  }
  gRNA.table[picked,]
}

pickedgRNAs.bygene <- list()

for (i in 1:length(uniquegRNAs.bygene)) {
  if (names(uniquegRNAs.bygene)[i] %in% names(selectgRNAs.bygene)) {
    pickedgRNAs.bygene[[i]] <- pickgrnas(selectgRNAs.bygene[[which(names(selectgRNAs.bygene) == names(uniquegRNAs.bygene)[i])]])
  }
}

names(pickedgRNAs.bygene) <- names(uniquegRNAs.bygene)

#Genes for which there are 5 or fewer unique gRNAs total are obviously going to ultimately choose all of the available gRNAs.
pickedgRNAs.bygene[which(sapply(uniquegRNAs.bygene, nrow) < 6)] <- uniquegRNAs.bygene[which(sapply(uniquegRNAs.bygene, nrow) < 6)]



#Why is YAL001C missing from selectgRNAs????
#You could probably slim down the filelist a lot here.

########
#TRASH
########

# #I then extracted the strain names with the following from the terminal:
# # head -n 59 1011Matrix.gvcf | tail -n 1 > ~/data/others_data/joseph_1011/strain.list
# #and then manually deleting the first 9 entries in textedit.
# #I then extracted the chromosome and position of each entry with the following in the terminal:
# # awk '{print $1,$2}' 1011Matrix.gvcf > snp.coords
# # tail +60 snp.coords > snp.coords.2
# 
# #Extracting the right part of the vcf (assuming there are any SNPs at all):
# #The code below takes:
# # (1) gstart.ind, corresponding to the index of the first variant position within the gRNA,
# # (2) gend.ind, corresponding to the index of the last variant position within the gRNA,
# # (3) Npos, the index of the N in the NGG relative to the index of the first variant position,
# #and returns which strains are perfect matches within that sequence, and which are not.
# 
# as.data.frame(
#   matrix(
#     scan("~/data/others_data/joseph_1011/1011Matrix.gvcf", 
#          what = "character", 
#          skip = 58 + gstart.ind, ###the index of the first snp in the range
#          nlines = 1 + gend.ind - gstart.ind, ###The number of snps in the range
#          ),
#     nrow = 1 + gend.ind - gstart.ind, ###The number of snps in the range
#     byrow = T)[,-c(1:9)],
#   stringsAsFactors = F) -> grna.vcf
# 
# if (Npos) grna.vcf[-Npos] -> grna.vcf #Remove the line corresponding to the N in the NGG
# 
# apply(grna.vcf, 2, function (x) min(sapply(x, function (y) substr(y, 1, 3) == "0/0"))) -> grna.perfect
#   
