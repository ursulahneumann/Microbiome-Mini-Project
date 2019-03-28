# Following DADA2 Pipeline Tutorial on a Small Multi-Sample Dataset
# Ursula H. Neumann
# Mar 25/19

# Introduction ----
# Tutorial instructions can be found here:
# https://benjjneb.github.io/dada2/tutorial.html

# Our starting point is a set of Illumina-sequenced paired-end fastq files that 
# have been split (or "demultiplexed") by sample and from which the 
# barcodes/adapters have already been removed. The end product is an amplicon 
# sequence variant (ASV) table, a higher-resolution analogue of the traditional 
# OTU table, which records the number of times each exact amplicon sequence 
# variant was observed in each sample. We also assign taxonomy to the output 
# sequences, and demonstrate how the data can be imported into the popular 
# phyloseq R package for the analysis of microbiome data.

# Set up ----
# Load packages
library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

# Set ggplot2 theme
theme_set(theme_bw())

# Define path and view files
path <- "C:/Users/ursula/Dropbox/Personal_Coding_Projects/Microbiome/dada2_tutorial/data"
list.files(path)

# Now we read in the names of the fastq files, and perform some string 
# manipulation to get matched lists of the forward and reverse fastq files
# String manipulation may need to be different for your own files

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and 
# SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Read quality ----
# Inspect quality of forward reads
# In gray-scale is a heat map of the frequency of each quality score at each 
# base position. The median quality score at each position is shown by the green
# line, and the quartiles of the quality score distribution by the orange lines. 
# The red line shows the scaled proportion of reads that extend to at least that
# position
plotQualityProfile(fnFs[1:2])

# Inspect quality of reverse reads.
plotQualityProfile(fnRs[1:2])

# Filter and trim ----
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Note - on Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)

# View first few lines of 'out'
head(out)

# Learn error rates ----
# Learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Visualize the estimated error rates
# The error rates for each possible transition (A-->C, A-->G, ...) are shown. 
# Points # are the observed error rates for each consensus quality score. The 
# black line shows the estimated error rates after convergence of the machine-
# learning algorithm. The red line shows the error rates expected under the 
# nominal definition of the Q-score. 
plotErrors(errF, nominalQ=TRUE)

# Dereplication ----
# Perform dereplication.  Dereplication combines all identical sequencing reads 
# into into "unique sequences" with a corresponding "abundance" equal to the 
# number of reads with that unique sequence. Dereplication substantially reduces
# computation time by eliminating redundant comparisons.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference ----
# Apply the core sample inference algorithm to the dereplicated data.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Inspect the returned dada-class object.
dadaFs[[1]]

# Merge paired reads ----
# Merge the forward/reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Sequence table ----
# Construct an amplicon sequence variant table (ASV) table.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ----
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through pipeline ----
# As a final check of our progress, we'll look at the number of reads that made 
# it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: 
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy ----
# Assign taxonomy to the sequence variants
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/ursula/Dropbox/Personal_Coding_Projects/Microbiome/dada2_tutorial/silva_nr_v128_train_set.fa.gz", multithread=TRUE)

# Extensions: The dada2 package also implements a method to make species level 
# assignments based on exact matching between ASVs and sequenced reference 
# strains. Recent analysis suggests that exact matching (or 100% identity) is 
# the only appropriate way to assign species to 16S gene fragments. Currently, 
# species-assignment training fastas are available for the Silva and RDP 16S 
# databases. 
taxa <- addSpecies(taxa, "C:/Users/ursula/Dropbox/Personal_Coding_Projects/Microbiome/dada2_tutorial/silva_species_assignment_v128.fa.gz")

# Inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Alternatives: The recently developed IdTaxa taxonomic classification method is 
# also available via the DECIPHER Bioconductor package. The paper introducing 
# the IDTAXA algorithm reports classification performance that is better than 
# the long-time standard set by the naive Bayesian classifier. 

# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim)) 
load("C:/Users/ursula/Dropbox/Personal_Coding_Projects/Microbiome/dada2_tutorial/silva_SSU_r132_March2018.RData") 
# processors=NULL or =2 causes my RStudio to crash
ids <- IdTaxa(dna, trainingSet, strand="top", processors=1, verbose=TRUE) 
# Ranks of interest
ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") 
# Convert the output object of class "Taxa" to a matrix analogous to the output 
# from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

# Redefine taxa to be taxid and redo analysis
taxid.print <- taxid # Removing sequence rownames for display only
rownames(taxid.print) <- NULL
head(taxid.print)

# Accuracy ---- 
# Evaluating DADA2's accuracy on the mock community
unqs.mock <- seqtab.nochim["Mock",]
# Drop ASVs absent in the Mock
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) 
cat("DADA2 inferred", length(unqs.mock), 
    "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), 
    "were exact matches to the expected reference sequences.\n")

# phyloseq plots ----
# Construct a simple sample data.frame from the information encoded in the 
# filenames
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

# Construct a phyloseq object directly from the dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

# Visualize alpha-diversity:
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

# Save Session ----
save.image("C:/Users/ursula/Dropbox/Personal_Coding_Projects/Microbiome/dada2_tutorial/dada2_tutorial.RData")