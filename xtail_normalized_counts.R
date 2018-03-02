rm(list = ls(all.names = TRUE))

library(xtail)

setwd("~/mcf7-ribo/analysis/results_R/")

# read ribosome footprint table with normalized read counts
rpf <- read.csv("norm_counts_cds_4_FP_samples.csv", row.names = 1)

# read total mRNA table with normalized read counts
mrna <- read.csv("norm_counts_cds_4_Total_samples.csv", row.names = 1)

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
write.csv(test_tab, "xtail_results_cds_4_samples.csv")

