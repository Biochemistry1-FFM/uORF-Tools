#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-b", "--bam_directory_path"), type = "character", default = NULL,
              help = "Path to directory containing .bam files", metavar = "character"),
  make_option(c("-a", "--uORF_annotation_file_path"), type = "character", default = NULL,
              help = "Path to .bed file with uORF annotation", metavar = "character"),
  make_option(c("-s", "--size_in_path"), type = "character", default = NULL,
              help = "Path for input size factor file", metavar = "character"),
  make_option(c("-t", "--sample_file_path"), type = "character", default = NULL,
              help = "Path to sample.tsv", metavar = "character"),
  make_option(c("-n", "--norm_count_uorf_out_path"), type = "character", default = "NULL",
              help = "Path for writing uORF normalized count file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$bam_directory_path)){
  print_help(option_parser)
  stop("Please supply arguments (-b, -a, -s, -t, -n), see --help \n", call.=FALSE)
}

library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)

# import uORFs bed file
uORFs <- import(options$uORF_annotation_file_path, format = "BED")

# create empty data frame, keep name as unique identifier
gene.counts <- data.frame(ORF.id = uORFs$name)

# get sample files
sample.files <- paste(options$bam_directory_path, list.files(options$bam_directory_path, pattern = "\\.bam$"), sep = "")

for (i in sample.files) {
  
  # get sample name
  name.i <- as.character(i)
  
  # import reads
  reads <- readGAlignments(i)
  
  # convert to granges
  reads <- granges(reads)
  
  # keep only first nt
  reads <- flank(reads, -1)
  
  # count reads into genes
  gene.counts[, name.i] <- countOverlaps(uORFs, reads)
  
}


# change row names and drop column gene.id
rownames(gene.counts) <- gene.counts$ORF.id
gene.counts$ORF.id <- NULL

# get sample sheet
sampleSheet <- read.csv(file=options$sample_file_path ,header=TRUE, sep="\t", stringsAsFactors=FALSE)

# col names and row names 
sampleName <- function(x) {
  method <- x[1]
  condition <- x[2]
  replicate <- x[3]
  sname <- paste(method, condition, replicate, sep="-")
  return(sname)
}

# generate sampleTable
sampleNames <- apply(sampleSheet,1,sampleName)
sample.files <- paste(sampleNames, ".bam", sep="")
conditions <- sampleSheet[,2]
sampleTable <- data.frame(row.names = sampleNames, fileName = sample.files, condition = conditions)

colnames(gene.counts) <- rownames(sampleTable)


# create DESeq data set
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = sampleTable,
                              design = ~ condition)


# supply size factors from whole library on longest protein coding
size.factors <- read.csv(options$size_in_path,row.names = 1, stringsAsFactors = FALSE)
colnames(size.factors) <- "size"

# apply size factors to dds object
sizeFactors(dds) <- size.factors$size

# get normalized counts
norm.counts <- counts(dds, normalized = TRUE)

# save normalized counts
write.csv(norm.counts, options$norm_count_uorf_out_path, quote = F)
