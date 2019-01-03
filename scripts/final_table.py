#!/usr/bin/env python
'''This script takes 2 input csv files generated by xtail_normalized_counts.R,
(one for uORFs and one for CDS) and ORF annotation from merged ribotish output
and computes the ribo_change parameter for regulatory uORFs and their mainORF.
'''

import pandas as pd
import argparse
import math


def set_change_symbol(log2change):
    if log2change > 1:
        return "+"
    else:
        return "-"


def uORF_change(uORFrowIn, ORFreadsIn):
    uORFrow = uORFrowIn
    ORFreads = ORFreadsIn
    replicates = math.ceil(len(uORFrow)/2)
    ratio1sum = 0
    ratio2sum = 0
    changesum = 0
    for replicate in range(0, replicates):
        uORFCond1 = uORFrow[replicate] + 1
        orfCond1 = ORFreads[replicate] + 1
        uORFCond2 = uORFrow[replicate + replicates] + 1
        orfCond2 = ORFreads[replicate + replicates] + 1
        ratio1 = orfCond1 / uORFCond1
        ratio2 = orfCond2 / uORFCond2
        change = ratio1 / ratio2
        ratio1sum += ratio1
        ratio2sum += ratio2
        changesum += change
    averageratio1 = ratio1sum / replicates 
    averageratio2 = ratio2sum / replicates
    averagechange = changesum / replicates
    return (averagechange,averageratio1,averageratio2)

def uORF_changes(uorf_table, uorf_reads_dict, orf_reads_dict):
    changes = []
    for _, uORFrow in uorf_table.iterrows():
        uORFid = uORFrow['uORFids']
        ORFid = uORFrow['transcript_id']
        uORFreads = uorf_reads_dict[uORFid]
        ORFreads = orf_reads_dict[ORFid]
        (averagechange,averageratio1,averageratio2) = uORF_change(uORFreads, ORFreads)
        joined_row = '\t'.join(map(str, uORFrow))
        logaveragechange = math.log2(averagechange)
        uORF_changes_string = joined_row + "\t" + str(averageratio1) + "\t" + str(averageratio2) + "\t" + str(averagechange) + "\t" + str(logaveragechange) + "\t" + set_change_symbol(averagechange)
        changes.append(uORF_changes_string)
    return (changes)


# read in xtail output files
def create_table(name):
    df = pd.read_table(name, sep=",", index_col=0)
    df = df[["log2FC_TE_final", "pvalue_final", "pvalue.adjust"]]
    return df


# create output data frame
def create_output(args):
    annot = pd.read_table(args.uORF_annotation, sep=",", index_col=0)
    annot.drop(columns=["chromosome", "start", "stop", "strand", "gene_id", "strand", "ORF_length"], axis=1, inplace=True)
    return annot


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Merges xtail differential analysis of translation efficiency of uORFs and their associated mainORF.')
    parser.add_argument("--uORF_reads", help='Path to input file with uORF reads')
    parser.add_argument("--ORF_reads", help='Path to input file with ORF reads')
    parser.add_argument("--uORF_annotation", help='Path to csv file containing uORF annotation')
    parser.add_argument("--output_csv_filepath", help='Path to write merged csv output')
    args = parser.parse_args()
    uorf_reads = pd.read_csv(args.uORF_reads)
    uorf_reads = uorf_reads[uorf_reads.columns.drop(list(uorf_reads.filter(regex='RNA')))]
    uorf_cols = uorf_reads.columns.values
    uorf_cols[0] = 'ID'
    uorf_reads.columns = uorf_cols
    uorf_reads_dict = uorf_reads.set_index('ID').T.to_dict('list')
    orf_reads = pd.read_csv(args.ORF_reads)
    orf_cols = orf_reads.columns.values
    orf_cols[0] = 'ID'
    orf_reads.columns = orf_cols
    orf_reads_dict = orf_reads.set_index('ID').T.to_dict('list')
    df_final = create_output(args)
    changes_list = uORF_changes(df_final, uorf_reads_dict, orf_reads_dict)
    changes_header = "coordinates\tgene_symbol\tstart_codon\ttranscript_id\tuORF_id\torf_uorf_ratio_c1\torf_uorf_ratio_c2\tribo_change\tlog_ribo_change\tregulation\n"
    changes_string = changes_header + '\n'.join(map(str, changes_list))
    f = open(args.output_csv_filepath, 'wt', encoding='utf-8')
    f.write(changes_string)


if __name__ == '__main__':
    main()
