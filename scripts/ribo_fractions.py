#!/usr/bin/env python
'''Computes differences between riboseq reads between conditions
'''

import pandas as pd
import re
import argparse
import numpy as np
import os

# function to create final data frame
def create_fractions(uORFs,ORFs):
    for index, row in df.uORFs():

    df_final="test"
    return df_final

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Computes differences between riboseq reads between conditions')
    parser.add_argument("--uORF_reads", help='Path to input file with uORF reads')
    parser.add_argument("--ORF_reads", help='Path to input file with ORF reads')
    parser.add_argument("--fraction_output", help='Path to write output file')
    args = parser.parse_args()
    uORFs = df = create_table(args.uORF_reads)
    ORFs = df = create_table(args.ORF_reads)
    ORFdictionary=ORFs.to_dict('list')
    fraction_output = create_fractions(uORFs,ORFs)
    f = open(args.fraction_output, 'wt', encoding='utf-8')
    f.write(fraction_output)
if __name__ == '__main__':
    main()
