#!/usr/bin/env python
'''Computes differences between riboseq reads between conditions
'''

import pandas as pd
import re
import argparse
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import math

def set_change_symbol(log2change):
    if log2change > 0: 
        return "+"
    else: 
        return "-"

def uORF_change(uORFrow,ORFreads): 
    print("uORFrow" + "\t"  + '\t'.join(map(str,uORFrow)) + "\n")
    print("ORFreads" + "\t" + '\t'.join(map(str,ORFreads))  + "\n")
    replicates = math.floor(len(uORFrow)/2)
    logchanges = []
    for replicate in range(1,replicates):
        uORFCond1 = uORFrow[replicate] + 1
        orfCond1 = ORFreads[replicate] + 1
        uORFCond2 = uORFrow[replicate + replicates] + 1
        orfCond2 = ORFreads[replicate + replicates] + 1
        ratio1 = orfCond1 / uORFCond1 
        ratio2 = orfCond2  / uORFCond2 #(WT, ref)
        change =  ratio1 / ratio2
        log2change = math.log2(change)
        logchanges.append(ratio1)
        logchanges.append(ratio2)
        logchanges.append(change)
        logchanges.append(set_change_symbol(log2change))
    logchangestring = '\t'.join(map(str,logchanges))
    print(logchangestring)
    return logchangestring
      
def uORF_changes(uORFreads,orfReadsDict):
    changes = []
    for _ , uORFrow in uORFreads.iterrows():
       uORFid = uORFrow[0]
       ORFid = re.sub(r'.\d+$', '', uORFid)
       ORFreads = orfReadsDict[ORFid]
       uORF_changes = uORF_change(uORFrow,ORFreads)
       changes.append(uORF_changes)
       #break
    return changes

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Computes differences between riboseq reads between conditions')
    parser.add_argument("--uORF_reads", help='Path to input file with uORF reads')
    parser.add_argument("--ORF_reads", help='Path to input file with ORF reads')
    parser.add_argument("--changes_output", help='Path to write output file')
    args = parser.parse_args()
    uORFreads = pd.read_csv(args.uORF_reads)
    ORFs = pd.read_csv(args.ORF_reads)
    new_columns = ORFs.columns.values
    new_columns[0] = 'ID'
    ORFs.columns = new_columns
    orfReadsDict=ORFs.set_index('ID').T.to_dict('list')
    changes_output = uORF_changes(uORFreads,orfReadsDict)
    changes_string = '\n'.join(map(str,changes_output))
    f = open(args.changes_output, 'wt', encoding='utf-8')
    f.write(changes_string)
if __name__ == '__main__':
    main()
