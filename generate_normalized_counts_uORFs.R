#!/usr/bin/env Rscript

library("optparse")

option_list = list(
  make_option(c("-b", "--bam_directory_path"), type = "character", default = NULL,
              help = "Path to directory containing .bam files", metavar = "character"),
  make_option(c("-s", "--size_in_path"), type = "character", default = NULL,
              help = "Path for input size factor file", metavar = "character"),
  make_option(c("-u", "--uorf_result_file_path"), type = "character", default = NULL,
              help = "Path to .csv file with merged uORF results", metavar = "character"),
  make_option(c("-n", "--norm_count_uorf_out_path"), type = "character", default = "NULL",
              help = "Path for writing uORF normalized count file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$bam_directory_path)){
  print_help(opt_parser)
  stop("Please supply arguments (-b, -a, -s, -n), see --help \n", call.=FALSE)
}


#rm(list = ls(all.names = TRUE))
#setwd('~/mcf7-ribo/data/seqs/bam_files_anica/')

library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(plyr)

# import uORFs table
#table <- read.csv("../../../analysis/ribotaper_anica/Merged_uORF_results.csv", header = TRUE, row.names = 1)
table <- read.csv(options$uorf_result_file_path, header = TRUE, row.names = 1)
table$gene_id <- as.character(table$gene_id)
table$transcript_id <- as.character(table$transcript_id)
table$gene_symbol <- as.character(table$gene_symbol)
table$ORF_id_gen <- as.character(table$ORF_id_gen)

# create uORFs ranges and adjust length column as these ranges include introns again
uORFs <- makeGRangesFromDataFrame(table, keep.extra.columns = TRUE)
uORFs$ORF_length <- width(uORFs)

# define sample type (RIBO ("FP_") or RNA ("Total_"))
sample.type <- "FP_"

# define bam file folder
#bam.folder <- '~/mcf7-ribo/data/seqs/bam_files_anica/'

# create empty data frame, keep ORF_id_gen as unique identifier (gene.ids in this case are not unique!)
gene.counts <- data.frame(ORF.id = uORFs$ORF_id_gen)

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
  gene.counts[, name.i] <- countOverlaps(uORFs, reads)

}

# change row names and drop column gene.id
rownames(gene.counts) <- gene.counts$ORF.id
gene.counts$ORF.id <- NULL

# set up sample table
#condition <- c(rep("control", each = 4), rep("treat", each = 4))
condition <- c(rep("control", each = 1), rep("treat", each = 1))
sampleTable <- data.frame(row.names = sample.names, fileName = sample.files,
                          condition = condition)
colnames(gene.counts) <- rownames(sampleTable)

# create DESeq data set
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = sampleTable,
                              design = ~ condition)


# supply size factors from whole library on longest protein coding
#size.factors <- read.csv("~/mcf7-ribo/analysis/results_R/size_factors_longest_protein_4_FP_samples.csv",row.names = 1, stringsAsFactors = FALSE)
size.factors <- read.csv("options$size_in_path",row.names = 1, stringsAsFactors = FALSE)
colnames(size.factors) <- "size"
sizeFactors(dds) <- size.factors$size

# differential analysis
dds <- DESeq(dds)

# get normalized counts
norm.counts <- counts(dds, normalized = TRUE)

# save normalized counts, change file name
#write.csv(norm.counts, "~/mcf7-ribo/analysis/results_R/norm_counts_uORFs_4_FP_samples.csv")
write.csv(norm.counts, options$norm_count_uorf_out_path)
