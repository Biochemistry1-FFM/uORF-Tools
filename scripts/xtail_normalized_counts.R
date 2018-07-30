#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-r", "--normalized_read_counts_csv_path"), type = "character", default = NULL,
              help = "Path to normalized read counts table", metavar = "character"),
  make_option(c("-t", "--sample_file_path"), type = "character", default = NULL,
              help = "Path to sample.tsv", metavar = "character"),
  make_option(c("-x", "--xtail_result_path"), type = "character", default = "NULL",
              help = "Path for writing xtail result file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$normalized_read_counts_csv_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -t, -x), see --help \n", call.=FALSE)
}

library(xtail)

# read table with mormalized read counts
counts <- read.csv(options$normalized_read_counts_csv_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# get sample sheet
sampleSheet <- read.csv(file=options$sample_file_path ,header=TRUE, sep="\t", stringsAsFactors=FALSE)

# split data frame into RIBO and RNA
RIBO <- counts[, sampleSheet$method == "RIBO"]
RNA <- counts[, sampleSheet$method == "RNA"]

#create condition vextor
condition <- sampleSheet$condition[which(sampleSheet$method == "RIBO")]

# run xtail analysis
test_results <- xtail(RNA, RIBO, condition, normalize = FALSE)

# turn results into table
test_tab <- resultsTable(test_results)

# write results into file
write.csv(test_tab, options$xtail_result_path, quote = F)

