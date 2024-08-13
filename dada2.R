library(dada2)
path <- "/core/cbc/fecal_metagenomics/scripts/amplicon_analysis/01_QC/trimmed_sequences" #fastq files,  V4 region of the 16S rRNA gene
#list.files(path) #view files 

#Forward and reverse fastq filenames have format: SAMPLENAME_trim_1.fastq and SAMPLENAME_trim_2.fastq
fnFs <- sort(list.files(path, pattern="_trim_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_trim_2.fastq", full.names = TRUE))
#Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#visualize the quality profiles of the forward reads
png(filename="QC1.png")
plotQualityProfile(fnFs[1:2])
dev.off()
#visualize quality profile of the reverse reads
png(filename="QC2.png")
plotQualityProfile(fnRs[1:2])
dev.off()

#place filtered files in filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#filter 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) 
head(out)

#error rates 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
png(filename="errate.png")
plotErrors(errF, nominalQ=TRUE)

#sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#Inspect the merger data.frame from the first sample
head(mergers[[1]])

#sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#accuracy
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

