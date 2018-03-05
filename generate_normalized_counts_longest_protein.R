#!/usr/bin/env Rscript

library("optparse")

option_list = list(
  make_option(c("-b", "--bam_directory_path"), type = "character", default = NULL,
              help = "Path to directory containing .bam files", metavar = "character"),
  make_option(c("-a", "--annotation_file_path"), type = "character", default = NULL,
              help = "Path to .gtf file with annotation", metavar = "character"),
  make_option(c("-s", "--size_out_path"), type = "character", default = NULL,
              help = "Path for writing output size file", metavar = "character"),
  make_option(c("-n", "--norm_count_out_path"), type = "character", default = "NULL",
              help = "Path for writing output normalized count file", metavar = "character")
);

if (is.null(opt$bam_directory_path)){
  print_help(opt_parser)
  stop("Please supply arguments (-b, -a, -s, -n), see --help \n", call.=FALSE)
}

option_parser = OptionParser(option_list=option_list);
options = parse_args(opt_parser);

#rm(list = ls(all.names = TRUE))
#setwd('~/mcf7-ribo/data/seqs/bam_files_anica/')
#use cwd instead

library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(plyr)

# import longest protein coding transcripts
#gencode <- import.gff("/shared/Homo_sapiens/NCBI/GRCh38/Annotation/Genes.gencode/gencode.v27.longest_protein_coding_transcript.gtf")
gencode <- import.gff(options$annotation_file_path)

# define sample type (RIBO ("FP_") or RNA ("Total_"))
sample.type <- "FP_"

# define bam file folder
#bam.folder <- '~/mcf7-ribo/data/seqs/bam_files_anica/'

# create empty data frame
gene.counts <- data.frame(gene.id = gencode$transcript_id)

# get sample files
sample.files <- paste(options$bam_directory_path, grep("FP_",list.files(options$bam_directory_path), value = TRUE), sep = "")

# exclue samples from experiment no. 1, keep re-sequencing experiment i.e. 1-2 (for now)
#sample.files <- sample.files[c(1,3,4,5,6,8,9,10)]
sample.files <- sample.files[c(1,2)]

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
  gene.counts[, name.i] <- countOverlaps(gencode, reads)

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

# differential analysis
dds <- DESeq(dds)
res <- results(dds)

# sizeFactors and normalized counts
size.factors <- sizeFactors(dds)
norm.counts <- counts(dds, normalized = TRUE)

# save normalized counts and size factors
#write.csv(size.factors, "~/mcf7-ribo/analysis/results_R/size_factors_logest_protein_4_FP_samples.csv")
#write.csv(norm.counts, "~/mcf7-ribo/analysis/results_R/norm_counts_longest_protein_4_FP_samples.csv")
write.csv(size.factors, "options$size_out_path")
write.csv(norm.counts, "options$norm_count_out_path")
