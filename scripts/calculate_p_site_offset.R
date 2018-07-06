#!/usr/bin/env Rscript

library(optparse)
library(plyr)

option_list = list(
  make_option(c("-i", "--metaplot_file_path"), type = "character", default = NULL,
              help = "Path to metaplot file", metavar = "character"),
  make_option(c("-o", "--offset_out_path"), type = "character", default = NULL,
              help = "Path for writing output offset file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser)


###script calculating P-site offset for ribotaper pipeline, takes as arguments ONE file created from create_metaplots.bash (table containing downsized bam-file info, in this case named RIBO-THAP-rep1)

# read in table created by create_metaplots.bash
reads<-read.table(options$metaplot_file_path,stringsAsFactors=F,header=F,sep="\t",comment.char="")

# change colnames
colnames(reads)<-c("chr","start","end","read_id","map_quality","strand",".1",".2",".3","spanning_exons","length_per_exon","length_introns","chr_stst","start_stst","end_stst","type_stst","gene_id_stst","strand_stst")

# select only reads without introns
reads_simpl<-reads[reads[,"length_introns"]=="0",]

# add column count with pseudocount of 1
reads_simpl$count<-1

# split data frame into list by strand
list_str<-split.data.frame(reads_simpl,f=reads_simpl[,"strand"])

# calculate distances
list_str[["+"]]$distance<-list_str[["+"]][,"start"]-list_str[["+"]][,"start_stst"]
list_str[["-"]]$distance<-list_str[["-"]][,"end_stst"]-list_str[["-"]][,"end"]

# merge lists into data frame
reads_simpl<-do.call(rbind.data.frame,list_str)

# count distance per read by length
dists_all<-with(reads_simpl,aggregate(count,by=list(type_stst,length_per_exon,distance),FUN=sum))

# change colnames
colnames(dists_all)<-c("type","length","distance","counts")

# select only start_codons 
starts <- dists_all[which(dists_all$type == "start_codon"),]

# get only entries upstream of start site
starts <- starts[which(starts$distance%in%c(min(starts$distance):0)),]

# split data fram by length
list_starts <- split.data.frame(starts, f=starts[,"length"])

# get length with max counts
length_max <- sort(sapply(list_starts, function(x) sum(x$counts)), decreasing = T)[c(1:5)]

# keep only 5 lengths with max counts
sel <- names(list_starts)%in%names(length_max)
list_max <- list_starts[which(sel)]

# get distance with max counts
distance <- sapply(list_max, function(x) x$distance[which(x$counts == max(x$counts))][1])

# create output data frame
result <- as.data.frame(t(distance))
result <- na.omit(result)

lengths <- names(result)
lengths <- na.omit(lengths)

lengthstring <- paste(as.character(lengths), collapse=",")
offsetstring <-paste(as.character(result), collapse=",")
output <- paste(lengthstring, offsetstring, sep=" ")

# write output
cat(output, file=options$offset_out_path, append=FALSE)

