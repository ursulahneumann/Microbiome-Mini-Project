# DADA2 Pipeline on Canine IBD samples 
# Ursula H. Neumann
# Mar 25/19

# Introduction ----
# Following dada2 tutorial found here:
# https://benjjneb.github.io/dada2/tutorial.html

# The canine samples only have forward reads so the analysis will be slightly
# different than what is shown in the tutorial which uses paired-end fastq 
# files.  I believe the barcodes and primers have already been removed.  

# I have been asked to investigate the healthy vs IBD only and this removed the
# acute hem. diarrhea samples (833.AHD.11 - 833.AHD.30).  

# The following samples had less than 15000 reads and were not included
# in the analysis.  This number may need to be changed. It may be 
# better to remove these samples programmatically in the future but
# I haven't figured that out yet.

# AHD.77
# AHD.79
# California.HC19
# California.HC23
# Finland.HC.16
# Finland.HC.18
# Finland.HC.5
# France.DC3
# Leda.A9
# Leda.A41
# Leda.A45
# Leda.A65
# Mel.TX.HC1
# Mel.TX.HC5
# Mel.TX.HC12
# Mel.TX.HC14
# Nor.C5
# Nor.C11
# Nor.C13
# Nor.C15
# RVC11
# RVC16
# Sweden.HC.18
# Sweden.IBD.063A
# Sweden.IBD.097A
# Sweden.IBD.103A

#           N (start)     N (reads>15000)
# total        177               150
# healthy      98                85
# IBD          79                65

# Set up ----
# Load packages
library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

# Set ggplot2 theme
theme_set(theme_bw())

# Define path and view files
path <- "C:/Users/ursula/Dropbox/Personal_Coding_Projects/Microbiome/dada2_canine_ibd/data"
list.files(path)

# Forward reads can be identified using the .fastq ending
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# Write function to extract sample names from filepaths so that it matches the
# anonymized_name column of the metadata file
extract_sample_ids <- function(filenames){
  processed <- unname(sapply(filenames, strsplit, '\\.')) # split on period
  processed <- sapply(processed, head, -1) # Remove the ".fastq"
  processed <- sapply(processed, tail, -1) # Remove the "833"
  processed <- sapply(processed, paste, collapse='.') # Paste with period
  return(processed)
}

# Extract sample names from filepaths
sample.names <- extract_sample_ids(fnFs)
sample.names

# Read quality ----
# Inspect quality of forward reads

# In gray-scale is a heat map of the frequency of each quality score at each 
# base position. The median quality score at each position is shown by the green
# line, and the quartiles of the quality score distribution by the orange lines. 
# The red line shows the scaled proportion of reads that extend to at least that
# position
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnFs[3:4])
plotQualityProfile(fnFs[5:6])

# Filter and trim ----
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

# Note - on Windows set multithread=FALSE
# A length of 75 was chosen based on the quality profile graphs. This should
# be investigated and tested further
out <- filterAndTrim(fnFs, filtFs, truncLen=c(75),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)

# View first few lines of 'out'
head(out)

# The code below was to try to detect and then remove any samples that had less
# than 1000 reads.  I ended up manually removing the files but would still be
# interested in trying to do this programmatically
#### Convert out to a df
#### out_df <- as.data.frame(out)
#### Determine which files have less than 1000 reads
#### files_low_reads <- row.names(out_df)[out_df$reads.in <1000]
#### Remove "833." from characters
#### files_low_reads <- substring(files_low_reads, first=5)
#### Generate path names to move files
#### files_low_reads <- paste(path, "/filtered/", files_low_reads, ".gz", sep="")

# Learn error rates ----
errF <- learnErrors(filtFs, multithread=TRUE)

# Visualize the estimated error rates
# The error rates for each possible transition (A-->C, A-->G, ...) are shown. 
# Points # are the observed error rates for each consensus quality score. The 
# black line shows the estimated error rates after convergence of the machine-
# learning algorithm. The red line shows the error rates expected under the 
# nominal definition of the Q-score. 
plotErrors(errF, nominalQ=TRUE)

# Dereplication ----
# Dereplication combines all identical sequencing reads into into "unique 
# sequences" with a corresponding "abundance" equal to the number of reads with 
# that unique sequence. Dereplication substantially reduces computation time by 
# eliminating redundant comparisons.
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

for (i in derepFs)
{
  print(min(i$quals))
}

# It looks like for many of the samples the minimum quality is 5 which seems 
# quite low.  The minimum values for both the forward and reverse reads in the
# tutorial were 12.  This may require further optimization.

# Sample inference ----
# Apply the core sample inference algorithm to the dereplicated data.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

# Inspect the returned dada-class object.
dadaFs[[1]]

# Sequence table ----
# Construct an amplicon sequence variant table (ASV) table.
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect the distribution of sequence lengths
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
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: 
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
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
# Below is from the tutorial which used a mock community.  I'm not sure this is
# relevant here.
#### Evaluating DADA2's accuracy on the mock community
##unqs.mock <- seqtab.nochim["Mock",]
#### Drop ASVs absent in the Mock
##unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) 
##cat("DADA2 inferred", length(unqs.mock), 
##    "sample sequences present in the Mock community.\n")

##mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
##match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
##cat("Of those,", sum(match.ref), 
##    "were exact matches to the expected reference sequences.\n")

# phyloseq plots ----

# Read in metadata file
metadata = read.table('C:/Users/ursula/Dropbox/Personal_Coding_Projects/Microbiome/dada2_canine_ibd/data/metadata.txt', 
                      sep='\t',
                      header = TRUE)

# Construct a simple sample data.frame from the information encoded in the 
# filenames and the metadata file
samples.out <- rownames(seqtab.nochim)
sample_boolean <- as.character(metadata$anonymized_name) %in% samples.out
metadata_incl <- metadata[sample_boolean,]

samdf <- data.frame(Disease_Status=metadata_incl$disease_stat)
rownames(samdf) <- metadata_incl$anonymized_name
samdf

# Construct a phyloseq object directly from the dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

# Visualize alpha-diversity
# This function gives a warning message suggesting to use un-trimmed data;
# however, the tutorial uses the function on the trimmed data.  I'm unsure
# whether we can ignore this error since we are using the dada2 method and not
# OTUs or if un-trimmed data should be used.
plot_richness(ps, measures=c("Shannon", "Simpson"), color="Disease_Status")
plot_richness(ps, measures=c("Chao1", "ACE"), color="Disease_Status")

# Graph alpha diversity with samples "grouped" by disease status
p <- plot_richness(ps, measures=c("Shannon", "Simpson"), color="Disease_Status")
new_order <- p$data$samples[order(p$data$Disease_Status[1:150])]
p$data$samples <- as.character(p$data$samples)
p$data$samples <- factor(p$data$samples, levels = new_order)
p

# Obtain raw measurements and merge with disease state column
alpha <- estimate_richness(ps, measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson","Fisher"))
alpha_df <- merge(alpha, samdf, by.x="row.names", by.y="row.names")

# Transform data to proportions as appropriate for Bray-Curtis distances
# This does not yet work as expected, will need to tweek.
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Disease_Status", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Family") + facet_wrap(~"Disease_Status", scales="free_x")

# Save Session ----
save.image("C:/Users/ursula/Dropbox/Personal_Coding_Projects/Microbiome/dada2_canine_ibd/dada2_canine_ibd.RData")