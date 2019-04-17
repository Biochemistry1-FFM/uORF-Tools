#!/usr/bin/env python
'''This script takes ORF annotation from merged ribotish output
and computes the ribo_change parameter for uORFs.
'''

import pandas as pd
import argparse
import math
import numpy as np

def create_replicateHeader(replicate_number):
    replicate_cond1 = []
    replicate_cond2 = []
    for replicate in range(0, replicate_number):
        cond1 = "Ratio-A-" + str(replicate)
        replicate_cond1.append(cond1)
        cond2 = "Ratio-B-" + str(replicate)
        replicate_cond2.append(cond2)
        cond1string = '\t'.join(map(str,replicate_cond1)) 
        cond2string = '\t'.join(map(str,replicate_cond2))
        replicate_header = cond1string + "\t" + cond2string
    return replicate_header

def uORF_change(uORFrowIn, ORFreadsIn):
    uORFrow = uORFrowIn
    ORFreads = ORFreadsIn
    replicate_number = math.ceil(len(uORFrow)/2)
    uorf1sum = 0
    orf1sum = 0
    uorf2sum = 0
    orf2sum = 0
    changesum = 0
    cond1_ratios = []
    cond2_ratios = []
    for replicate in range(0, replicate_number):
        uORFCond1 = uORFrow[replicate] + 1
        orfCond1 = ORFreads[replicate] + 1
        uORFCond2 = uORFrow[replicate + replicates] + 1
        orfCond2 = ORFreads[replicate + replicates] + 1
        ratio1 = orfCond1 / uORFCond1
        cond1_ratios.append(ratio1)
        ratio2 = orfCond2 / uORFCond2
        cond2_ratios.append(ratio2)
        change = ratio1 / ratio2
        uorf1sum += uORFCond1
        orf1sum += orfCond1
        uorf2sum += uORFCond2
        orf2sum += orfCond2
        changesum += change
    averageuORF1 = uorf1sum / replicates 
    averageORF1 = orf1sum / replicates
    averageuORF2 = uorf2sum / replicates 
    averageORF2 = orf2sum / replicates
    averagechange = changesum / replicates
    logaveragechange = math.log2(averagechange)
    return (cond1ratios,cond2ratios,logaveragechange)

def uORF_changes(uorf_table, uorf_reads_dict, orf_reads_dict):
    output = []
    for _, uORFrow in uorf_table.iterrows():
        uORFid = uORFrow['uORFids']
        ORFid = uORFrow['transcript_id']
        uORFreads = uorf_reads_dict[uORFid]
        ORFreads = orf_reads_dict[ORFid]
        (cond1ratios,cond2ratios,logaveragechange) = uORF_change(uORFreads, ORFreads)
        annotation_row = '\t'.join(map(str, uORFrow))
        cond1ratiosstring = '\t'.join(map(str,cond1ratios)) 
        cond2ratiosstring = '\t'.join(map(str,cond2ratios)) 
        uORF_changes_string=annotation_row + "\t" + cond1ratiosstring + "\t" + cond2ratiosstring + "\t" + logaveragechange
        output.append(uORF_changes_string)
    return (output)


# create output data frame
def create_output(args):
    annot = pd.read_csv(args.uORF_annotation, sep=",", index_col=0)
    annot.drop(columns=["chromosome", "start", "stop", "strand", "start_codon", "gene_id", "strand", "ORF_length"], axis=1, inplace=True)
    return annot


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Takes ORF annotation from merged ribotish output and computes the ribo_change parameter.')
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
    replicate_number = math.ceil(len(df_final[0])/2)
    changes_list = uORF_changes(df_final, uorf_reads_dict, orf_reads_dict)
    ratios_header = create_ratioheader(replicate_number)
    changes_header = "coordinates\tgene_symbol\ttranscript_id\tuORF_id\t" + ratios_header + "\tlog2FC_main_ORF_to_uORF_ratios\n"
    changes_string = changes_header + '\n'.join(map(str, changes_list))
    f = open(args.output_csv_filepath, 'wt', encoding='utf-8')
    f.write(changes_string)


if __name__ == '__main__':
    main()
