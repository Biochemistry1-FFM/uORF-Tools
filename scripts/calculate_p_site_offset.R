#!/usr/bin/Rscript

###script calculating P-site offset for ribotaper pipeline, takes as arguments ONE file created from create_metaplots.bash (table containing downsized bam-file info, in this case named RIBO-THAP-rep1)

print(paste("--- calculating P-site offset","---",date(),sep=" "))

args <- commandArgs(trailingOnly = TRUE)

# hardcoded file path, needs to be fixed
args <- c("../RIBO-THAP-rep1")

# read in table created by create_metaplots.bash
reads<-read.table(args[1],stringsAsFactors=F,header=F,sep="\t",comment.char="")

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

# select only start_codons and read lengths between 25-30
starts <- dists_all[which(dists_all$type == "start_codon" & dists_all$length%in%c(25:30)),]

# get only entries upstream of start site
starts <- starts[which(starts$distance%in%c(-23:0)),]

# split data fram by length
list_starts <- split.data.frame(starts, f=starts[,"length"])

# get distance with max counts
distance <- sapply(list_starts, function(x) x$distance[which(x$counts == max(x$counts))])

# create output data frame
output <- as.data.frame(t(distance))
output <- abs(output)

# write output
write.csv(output, "../offsets", row.names = F, quote = F)


