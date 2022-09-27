### DADA2 - 16 S analysis from paired-end illumina reads ###
# Script by: David Benito-Merino (dbenito@mpi-bremen.de).
# Based on the DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html

# DADA2 works on paired-end Illumina reads in which non-biological sequences (primers, adapters, etc.) have been removed.

require(dada2); require(Rcpp); require(ggplot2); require(phyloseq)

rm(list=ls())
dir <- "" # path to the directory containing the clipped files.
tmp <- ""
output <- ""
setwd(dir)
getwd()

load(paste(tmp, "dada2_archaea.RData", sep="/"))

list.files(dir)

# Forward and reverse fastq filenames have format: *_R1*.fastq and *_R2*.fastq
R1_list <- sort(list.files(dir, pattern="_R1", full.names = TRUE))
R2_list <- sort(list.files(dir, pattern="_R2", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.
# In our case the samples are named "1", "2", etc. Conserve the fileID_R1 and fileID_R2 to refer to original sample names.
sample.names_R1 <- sapply(strsplit(basename(R1_list), ".fastq"), `[`, 1)
sample.names_R2 <- sapply(strsplit(basename(R2_list), ".fastq"), `[`, 1)

# Visualisation of the reads' quality:
quality_plot_R1 <- plotQualityProfile(R1_list[1:length(R1_list)])
quality_plot_R2 <- plotQualityProfile(R2_list[1:length(R2_list)])

# Save plots
ggsave(filename="quality_plots_R1.pdf", path=output, plot=quality_plot_R1, device="pdf", height=30, width=30)
ggsave(filename="quality_plots_R2.pdf", path=output, plot=quality_plot_R2, device="pdf", height=30, width=30)

# Place filtered files in filtered/ subdirectory
filt_R1 <- file.path(path=dir, "filtered", paste0(sample.names_R1, "_R1_filt.fastq.gz"))
filt_R2 <- file.path(path=dir, "filtered", paste0(sample.names_R2, "_R2_filt.fastq.gz"))
names(filt_R1) <- sample.names_R1
names(filt_R2) <- sample.names_R2

# Standard filtering paremeters:
# Pick filtering parameters based on quality plots.
out <- filterAndTrim(R1_list, filt_R1, R2_list, filt_R2, truncLen=c(250,200),
              maxN=0, maxLen=300, minLen=100, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=T) # On Windows set multithread=FALSE
head(out)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))


# Calculate & plot error rates for nucleotide substitution
error_R1 <- learnErrors(filt_R1, multithread=T)
error_R2 <- learnErrors(filt_R2, multithread=T)

error_plot_R1 <- plotErrors(error_R1, nominalQ=T)
error_plot_R2 <- plotErrors(error_R2, nominalQ=T)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))


# Dereplication
derep_R1 <- derepFastq(filt_R1, verbose=T)
derep_R2 <- derepFastq(filt_R2, verbose=T)
# Use the same sample names:
names(derep_R1) <- sample.names_R1
names(derep_R2) <- sample.names_R2

# Sample inference algorithm (dada)
dada_R1 <- dada(derep_R1, err=error_R1, multithread=T, pool=T)
dada_R2 <- dada(derep_R2, err=error_R2, multithread=T, pool=T)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))

# Merging. Not supported for archaeal amplicons because the overlap is short.
# Non-overlapping read pairs can be merged with the parameter "justConcatenate=T" which adds a string of 10 N. This does not allow to assign species in merged reads.
# The solution is to assign species to forward and/or reverse reads separately.
nochim_R1 <- removeBimeraDenovo(dada_R1)
nochim_R2 <- removeBimeraDenovo(dada_R2)

# Get uniques for R1:
for(sam in sample.names_R1) {
  cat("Uniques:", sam, "\n")
  uniq_R1 <- nochim_R1
  uniq_R1[[sam]] <- getUniques(uniq_R1[[sam]])
}

# Get uniques for R2:
for(sam in sample.names_R2) {
  cat("Uniques:", sam, "\n")
  uniq_R2 <- nochim_R2
  uniq_R2[[sam]] <- getUniques(uniq_R2[[sam]])
}


# Amplicon sequence variants tables (ASV)
ASV_table_R1 <- makeSequenceTable(uniq_R1)
ASV_table_R2 <- makeSequenceTable(uniq_R2)

dim(ASV_table_R1)
dim(ASV_table_R2)


# Inspect distribution of sequence lengths:
seq_table <- table(nchar(getSequences(ASV_table)))
seq_hist <- hist(nchar(getSequences(ASV_table)))

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))

# Read counts throughout pipeline.
getN <- function(x) sum(getUniques(x))

track <- cbind(out,           +
               sapply(dada_R1, getN), +
               sapply(dada_R2, getN), +
               rowSums(ASV_table_R1), +
               rowSums(ASV_table_R2)
               )

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "nonchimF", "nonchimR")
rownames(track) <- sample.names_R1
head(track)

write.table(track, file=paste(output, "dada_read_track.tsv", sep="/"), sep="\t", quote=F)

# Taxonomic classification on forward and reverse reads separately:
# For R1: this works only with "assignSpecies" and not with "addSpecies"!
#for(sam in sample.names_R1) {
#  cat("Assign taxonomy:", sam, "\n")
#  # The database in the correct format can be found in the dada2 website.
#  ASV_taxonomy <- uniq_R1
#  ASV_taxonomy[[sam]] <- assignTaxonomy(uniq_R1[[sam]], "~/db/dada2/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=T)
#  cat("Add species:", sam, "\n")
#  # Add species:
#  ASV_taxonomy[[sam]] <- assignSpecies(uniq_R1[[sam]], "~/db/dada2/tax/silva_species_assignment_v138.1.fa.gz")
#  cat("Done", "\n", "\n")
#}

ASV_taxonomy_R1 <- assignTaxonomy(ASV_table_R1, "~/db/dada2/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=T)
ASV_taxonomy_R1 <- addSpecies(ASV_taxonomy_R1, "~/db/dada2/tax/silva_species_assignment_v138.1.fa.gz")
save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))


# For R2:
# The reads are in reverse complementary, so use the package "ShortRead" to transform them.
#for(sam in sample.names_R2) {
#  cat("Assign taxonomy:", sam, "\n")
#  # The database in the correct format can be found in the dada2 website.
#  ASV_taxonomy_R2 <- uniq_R2
#  ASV_taxonomy_R2[[sam]] <- assignTaxonomy(uniq_R2[[sam]], "~/db/dada2/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=T)
#  cat("Add species:", sam, "\n")
#  # Add species:
#  ASV_taxonomy_R2[[sam]] <- assignSpecies(uniq_R2[[sam]], "~/db/dada2/tax/silva_species_assignment_v138.1.fa.gz")
#  cat("Done", "\n", "\n")
#}

# Add the parameter tryRC=TRUE for reverse complementary reads (reverse read files are in reverse complementary)
ASV_taxonomy_R2 <- assignTaxonomy(ASV_table_R2, "~/db/dada2/tax/silva_nr99_v138.1_train_set.fa.gz", tryRC = TRUE, multithread=T)
ASV_taxonomy_R2 <- addSpecies(ASV_taxonomy_R2, "~/db/dada2/tax/silva_species_assignment_v138.1.fa.gz")
save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))

########################
# Analyse and plot results with phyloseq.
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

### For R1:
# Extract the sample and ASV names:
samples.out <- rownames(ASV_table_R1)
ASVs_R1 <- colnames(ASV_table_R1)
# ASVs ID table:
ASVs_ID_R1 <- cbind(ASVs_R1, paste("asv", c(1:ncol(ASV_table_R1)), sep=""))


# rename the ASV to asv#:
colnames(ASV_table_R1) <- paste("asv", c(1:ncol(ASV_table_R1)), sep="")
rownames(ASV_taxonomy_R1) <- paste("asv", c(1:nrow(ASV_taxonomy_R1)), sep="")
ASV_taxonomy_R1[is.na(ASV_taxonomy_R1[,1])] <- "Unclassified" # Replace empty taxons (domain/kingdom level) with "Unclassified".


# Add sample names:
head (samples.out)
samples.out <- cbind(samples.out, c("Hexadecane37", "Control37", "Hexadecane70",
                                    "Control70_old", "Hexadecane50", "Control70",
                                    "Sediment_4869_subset", "Sediment_4869", "Sediment_4869_0-2_subset", "Sediment_4869_0-2 ",
                                    "Sediment_4869_6-10_subset", "Sediment_4869_6-10", "Sediment_combined")
                    )
colnames(samples.out) <- c("ID", "Sample")
rownames(samples.out) <- samples.out[,1] # Row names are samples IDs.
samples.out <- as.data.frame(samples.out)
# replace in ASV_nochim:
rownames(ASV_table_R1) <- samples.out[,2]
rownames(samples.out) <- samples.out[,2]

# Create otu table:
OTU_phyloseq_R1 <- otu_table(ASV_table_R1, taxa_are_rows = FALSE)
SAMPLE_phyloseq_R1 <- sample_data(samples.out)
TAX_phyloseq_R1 <- tax_table(ASV_taxonomy_R1)

# Create a phyloseq object.
dada2_phyloseq_R1 <- phyloseq(OTU_phyloseq_R1, TAX_phyloseq_R1, SAMPLE_phyloseq_R1)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))

# Remove bacterial reads:
keep <- taxa_names(tax_table(dada2_phyloseq_R1)[tax_table(dada2_phyloseq_R1)[,1]=="Archaea"])
dada2_phyloseq_prune_R1 <- prune_taxa(keep, dada2_phyloseq_R1)

# Remove OTUs that occur less than 5 times in the whole dataset
dada2_phyloseq_prune_R1 <- prune_taxa(taxa_sums(dada2_phyloseq_prune_R1) > 5, dada2_phyloseq_prune_R1)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))


# Agglomerate taxa at order level:
dada2_phyloseq_prune_R1 <- tax_glom(dada2_phyloseq_prune_R1, "Order")

# To classify the NAs as a higher known taxon:
NameTax <- function(x, ind){
  if(is.na(x[ind])){
    x[ind] <- x[ind]
  } else {
    if(ind==1){x[ind] <- paste("k", x[ind], sep="__")} else{                  # Domain/Kingdom
      if(ind==2){x[ind] <- paste("p", x[ind], sep="__")} else{                # Phylum
        if(ind==3){x[ind] <- paste("c", x[ind], sep="__")} else{              # Class
          if(ind==4){x[ind] <- paste("o", x[ind], sep="__")} else{            # Order
            if(ind==5){x[ind] <- paste("f", x[ind], sep="__")} else{          # Family
              if(ind==6){x[ind] <- paste("g", x[ind], sep="__")} else{        # Genus
                if(ind==7){x[ind] <- paste("s", x[ind-1], x[ind], sep="__")}  # Species
              }
            }
          }
        }
      }
    }
  }
}


ModifyTax <- function(x,ind){
  #   xth row in the dataframe
  #   ind taxonomy level to change
  if(is.na(x[ind])){
    nonNa <- which(!is.na(x[-ind])) # which taxa are not NA excepting the one we're interested in.
    maxNonNa <- max(nonNa)
    x[ind] <- x[maxNonNa]
  }else{x[ind] <- x[ind]}
}

# Apply the funcitons NameTax and ModifyTax
tax.tab_R1 <- data.frame(tax_table(dada2_phyloseq_prune_R1))

for (i in 1:7) {
  tax_table(dada2_phyloseq_prune_R1)[,i] <- apply(tax.tab_R1, 1, NameTax, ind=i)
}

tax.tab_R1 <- data.frame(tax_table(dada2_phyloseq_prune_R1))

for (i in 1:7) {
  tax_table(dada2_phyloseq_prune_R1)[,i] <- apply(tax.tab_R1,1,ModifyTax,ind=i)
}

# Transform to relative abundance
physeq_RA_R1 <- transform_sample_counts(dada2_phyloseq_prune_R1, function(x) x/sum(x))

# Convert to dataframe:
data_R1 <- psmelt(physeq_RA_R1) # create dataframe from phyloseq object


# Simple way to rename  with < 1% abundance. Check dimensions of "data", i are the columns from phylum to species.
for (i in 7:12) {
data_R1[data_R1$Abundance < 0.01, i] <- "< 1% rel. abund."
}

# Order the samples
data_R1$Sample <- as.factor(data_R1$Sample)
data_R1$Sample <- factor(data_R1$Sample, levels=c("Sediment_4869_0-2 ", "Sediment_4869_0-2_subset", "Sediment_4869_6-10", "Sediment_4869_6-10_subset",
                             "Sediment_4869", "Sediment_4869_subset", "Sediment_combined", "Control37", "Control70_old", "Control70",
                             "Hexadecane37", "Hexadecane50",  "Hexadecane70"))

# Plots
plot_phylum_R1 <- ggplot(data_R1[order(data_R1$Phylum, decreasing=F), ], aes(x=Sample, y=Abundance, fill=factor(Phylum, levels=unique(Phylum))))+
                      geom_bar(aes(), stat="identity", position="stack", color="black")+
                      theme_classic()+
                      theme(axis.text.x = element_text(angle = 90))
plot_order_R1 <- ggplot(data_R1[order(data_R1$Phylum, decreasing=F), ], aes(x=Sample, y=Abundance, fill=factor(Order, levels=unique(Order))))+
                      geom_bar(aes(), stat="identity", position="stack", color="black")+
                      theme_classic()+
                      theme(axis.text.x = element_text(angle = 90))

# Save the plots:
ggsave(filename="plot_order_RA_R1.pdf", path=output, plot=plot_order_R1, device="pdf", height=9, width=16)
ggsave(filename="plot_phylum_RA_R1.pdf", path=output, plot=plot_phylum_R1, device="pdf", height=9, width=16)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))


#########

### For R2:

# Extract the sample and ASV names:
samples.out <- rownames(ASV_table_R2)
ASVs_R2 <- colnames(ASV_table_R2)
# ASVs ID table:
ASVs_ID_R2 <- cbind(ASVs_R2, paste("asv", c(1:ncol(ASV_table_R2)), sep=""))


# rename the ASV to asv#:
colnames(ASV_table_R2) <- paste("asv", c(1:ncol(ASV_table_R2)), sep="")
rownames(ASV_taxonomy_R2) <- paste("asv", c(1:nrow(ASV_taxonomy_R2)), sep="")
ASV_taxonomy_R2[is.na(ASV_taxonomy_R2[,1])] <- "Unclassified" # Replace empty taxons (domain/kingdom level) with "Unclassified".


# Add sample names:
head (samples.out)
samples.out <- cbind(samples.out, c("Hexadecane37", "Control37", "Hexadecane70",
                                    "Control70_old", "Hexadecane50", "Control70",
                                    "Sediment_4869_subset", "Sediment_4869", "Sediment_4869_0-2_subset", "Sediment_4869_0-2 ",
                                    "Sediment_4869_6-10_subset", "Sediment_4869_6-10", "Sediment_combined")
                    )
colnames(samples.out) <- c("ID", "Sample")
rownames(samples.out) <- samples.out[,1] # Row names are samples IDs.
samples.out <- as.data.frame(samples.out)
# replace in ASV_nochim:
rownames(ASV_table_R2) <- samples.out[,2]
rownames(samples.out) <- samples.out[,2]

# Create otu table:
OTU_phyloseq_R2 <- otu_table(ASV_table_R2, taxa_are_rows = FALSE)
SAMPLE_phyloseq_R2 <- sample_data(samples.out)
TAX_phyloseq_R2 <- tax_table(ASV_taxonomy_R2)

# Create a phyloseq object.
dada2_phyloseq_R2 <- phyloseq(OTU_phyloseq_R2, TAX_phyloseq_R2, SAMPLE_phyloseq_R2)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))

# Remove bacterial reads:
keep <- taxa_names(tax_table(dada2_phyloseq_R2)[tax_table(dada2_phyloseq_R2)[,1]=="Archaea"])
dada2_phyloseq_prune_R2 <- prune_taxa(keep, dada2_phyloseq_R2)

# Remove OTUs that occur less than 5 times in the whole dataset
dada2_phyloseq_prune_R2 <- prune_taxa(taxa_sums(dada2_phyloseq_prune_R2) > 5, dada2_phyloseq_prune_R2)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))


# Agglomerate taxa at order level:
dada2_phyloseq_prune_R2 <- tax_glom(dada2_phyloseq_prune_R2, "Order")

# To classify the NAs as a higher known taxon:
NameTax <- function(x, ind){
  if(is.na(x[ind])){
    x[ind] <- x[ind]
  } else {
    if(ind==1){x[ind] <- paste("k", x[ind], sep="__")} else{                  # Domain/Kingdom
      if(ind==2){x[ind] <- paste("p", x[ind], sep="__")} else{                # Phylum
        if(ind==3){x[ind] <- paste("c", x[ind], sep="__")} else{              # Class
          if(ind==4){x[ind] <- paste("o", x[ind], sep="__")} else{            # Order
            if(ind==5){x[ind] <- paste("f", x[ind], sep="__")} else{          # Family
              if(ind==6){x[ind] <- paste("g", x[ind], sep="__")} else{        # Genus
                if(ind==7){x[ind] <- paste("s", x[ind-1], x[ind], sep="__")}  # Species
              }
            }
          }
        }
      }
    }
  }
}


ModifyTax <- function(x,ind){
  #   xth row in the dataframe
  #   ind taxonomy level to change
  if(is.na(x[ind])){
    nonNa <- which(!is.na(x[-ind])) # which taxa are not NA excepting the one we're interested in.
    maxNonNa <- max(nonNa)
    x[ind] <- x[maxNonNa]
  }else{x[ind] <- x[ind]}
}

# Apply the funcitons NameTax and ModifyTax
tax.tab_R2 <- data.frame(tax_table(dada2_phyloseq_prune_R2))

for (i in 1:7) {
  tax_table(dada2_phyloseq_prune_R2)[,i] <- apply(tax.tab_R2, 1, NameTax, ind=i)
}

tax.tab_R2 <- data.frame(tax_table(dada2_phyloseq_prune_R2))

for (i in 1:7) {
  tax_table(dada2_phyloseq_prune_R2)[,i] <- apply(tax.tab_R2,1,ModifyTax,ind=i)
}

# Transform to relative abundance
physeq_RA_R2 <- transform_sample_counts(dada2_phyloseq_prune_R2, function(x) x/sum(x))

# Convert to dataframe:
data_R2 <- psmelt(physeq_RA_R2) # create dataframe from phyloseq object


# Simple way to rename  with < 1% abundance. Check dimensions of "data", i are the columns from phylum to species.
for (i in 7:12) {
data_R2[data_R2$Abundance < 0.01, i] <- "< 1% rel. abund."
}

# Order the samples
data_R2$Sample <- as.factor(data_R2$Sample)
data_R2$Sample <- factor(data_R2$Sample, levels=c("Sediment_4869_0-2 ", "Sediment_4869_0-2_subset", "Sediment_4869_6-10", "Sediment_4869_6-10_subset",
                             "Sediment_4869", "Sediment_4869_subset", "Sediment_combined", "Control37", "Control70_old", "Control70",
                             "Hexadecane37", "Hexadecane50",  "Hexadecane70"))

# Plots
plot_phylum_R2 <- ggplot(data_R2[order(data_R2$Phylum, decreasing=F), ], aes(x=Sample, y=Abundance, fill=factor(Phylum, levels=unique(Phylum))))+
                      geom_bar(aes(), stat="identity", position="stack", color="black")+
                      theme_classic()+
                      theme(axis.text.x = element_text(angle = 90))
plot_order_R2 <- ggplot(data_R2[order(data_R2$Phylum, decreasing=F), ], aes(x=Sample, y=Abundance, fill=factor(Order, levels=unique(Order))))+
                      geom_bar(aes(), stat="identity", position="stack", color="black")+
                      theme_classic()+
                      theme(axis.text.x = element_text(angle = 90))

# Save the plots:
ggsave(filename="plot_order_RA_R2.pdf", path=output, plot=plot_order_R2, device="pdf", height=9, width=16)
ggsave(filename="plot_phylum_RA_R2.pdf", path=output, plot=plot_phylum_R2, device="pdf", height=9, width=16)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))

#############
# final plot:

data_subset <- subset(data_R1, Sample %in% c("Sediment_combined", "Control37", "Control70", "Hexadecane37", "Hexadecane50", "Hexadecane70"))

plot_order_R1_final  <- ggplot(data_subset[order(data_subset$Phylum, decreasing=F), ],
                      aes(x=Sample, y=Abundance, fill=factor(Order, levels=unique(Order))))+
                      geom_bar(aes(), stat="identity", position="stack", color="black")+
                      theme_classic()+
                      guides(fill=guide_legend(ncol=2))+ # Split legend in two columns
                      theme(axis.text.x = element_text(angle = 90))

plot_order_R1_final_no_legend <- plot_order_R1_final + theme(legend.position = "none")

ggsave(filename="plot_order_RA_R1_final.pdf", path=output, plot=plot_order_R1_final_no_legend, device="pdf", height=9, width=6)
ggsave(filename="plot_order_RA_R1_final_legend.pdf", path=output, plot=plot_order_R1_final, device="pdf", height=9, width=6)

save.image(file=paste(tmp, "dada2_archaea.RData", sep="/"))
