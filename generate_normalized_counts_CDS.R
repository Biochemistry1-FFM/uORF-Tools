#!/usr/bin/env Rscript

rm(list = ls(all.names = TRUE))
setwd('~/mcf7-ribo/data/seqs/bam_files_anica/')

library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(plyr)

# import longest protein coding transcripts
gencode <- import.gff("/shared/Homo_sapiens/NCBI/GRCh38/Annotation/Genes.gencode/gencode.v27.longest_protein_coding_transcript.gtf")

# keep only cds
sel <- gencode$type == "CDS"
cds <- gencode[which(sel),]

# define sample type (RIBO ("FP_") or RNA ("Total_"))
sample.type <- "FP_"

# define bam file folder
bam.folder <- '~/mcf7-ribo/data/seqs/bam_files_anica/'

# create empty data frame
gene.counts <- data.frame(gene.id = cds$transcript_id)

# get sample files
sample.files <- paste(bam.folder, grep("FP_",list.files(bam.folder), value = TRUE), sep = "")

# exclue samples from experiment no. 1, keep re-sequencing experiment i.e. 1-2 (for now)
sample.files <- sample.files[c(1,3,4,5,6,8,9,10)]

# extract sample names
sample.names <- regmatches(sample.files,regexpr("FP_.*_[0-9]",sample.files))

for (i in sample.files) {
  
  # get sample name
  name.i <- regmatches(i,regexpr("FP_.*_[0-9]",i))
  
  # import reads
  reads <- readGAlignments(i)
  
  # get read lengths
  widths <- qwidth(reads)
  
  # convert to granges
  reads <- granges(reads)
  mcols(reads)$qwidth <- widths
  
  # keep only first nt
  reads <- flank(reads, -1)
  
  # keep only reads of 25-35 nt
  reads <- reads[elementMetadata(reads)$qwidth%in%c(25:35)]
  
  # count reads into genes
  gene.counts[, name.i] <- countOverlaps(cds, reads)
  
}

# sum all reads across gene.ids 
gene.counts <- ddply(gene.counts,"gene.id",numcolwise(sum))

# change row names and drop column gene.id
rownames(gene.counts) <- gene.counts$gene.id
gene.counts$gene.id <- NULL

# set up sample table
condition <- c(rep("control", each = 4), rep("treat", each = 4))
sampleTable <- data.frame(row.names = sample.names, fileName = sample.files,
                          condition = condition)
colnames(gene.counts) <- rownames(sampleTable)

# create DESeq data set
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = sampleTable,
                              design = ~ condition)


# supply size factors from whole library on longest protein coding
size.factors <- read.csv("~/mcf7-ribo/analysis/results_R/size_factors_logest_protein_4_FP_samples.csv",row.names = 1, stringsAsFactors = FALSE)
colnames(size.factors) <- "size"
sizeFactors(dds) <- size.factors$size

# differential analysis
dds <- DESeq(dds)

# get normalized counts
norm.counts <- counts(dds, normalized = TRUE)

# save normalized counts, change file name
write.csv(norm.counts, "~/mcf7-ribo/analysis/results_R/norm_counts_cds_4_FP_samples.csv")

