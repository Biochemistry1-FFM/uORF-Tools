#!/usr/bin/env Rscript

library(optparse)
library(xtail)
#rm(list = ls(all.names = TRUE))
option_list = list(
  make_option(c("-r", "--ribosome_footprint_normalized_read_counts_path"), type = "character", default = NULL,
              help = "Path to ribosome footprint table with normalized read counts", metavar = "character"),
  make_option(c("-n", "--normalized_read_counts_path"), type = "character", default = NULL,
              help = "Path to total mRNA table with normalized read counts", metavar = "character"),
  make_option(c("-x", "--xtail_result_path"), type = "character", default = "NULL",
              help = "Path for writing xtail result file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$ribosome_footprint_normalized_read_counts_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -n, -x), see --help \n", call.=FALSE)
}

#setwd("~/mcf7-ribo/analysis/results_R/")

# read ribosome footprint table with normalized read counts
#rpf <- read.csv("norm_counts_cds_4_FP_samples.csv", row.names = 1)
rpf <- read.csv(options$ribosome_footprint_normalized_read_counts_path, row.names = 1)

# read total mRNA table with normalized read counts
#mrna <- read.csv("norm_counts_cds_4_Total_samples.csv", row.names = 1)
mrna <- read.csv(options$normalized_read_counts_path, row.names = 1)


# create condition vector
condition <- c(rep("control", each = 4), rep("treat", each = 4))

# run xtail analysis
test_results <- xtail(mrna, rpf, condition, normalize = FALSE)

# turn reuslts into table
test_tab <- resultsTable(test_results)
test_tab <- test_tab[order(test_tab$pvalue.adjust),]

# some plots
plotFCs(test_results)
plotRs(test_results)
volcanoPlot(test_results)

# write results into file
#write.csv(test_tab, "xtail_results_cds_4_samples.csv")
write.csv(test_tab, options$xtail_result_path)
