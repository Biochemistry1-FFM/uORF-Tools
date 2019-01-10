#!/usr/bin/env python
'''Converts uORF annotation.csv to bed format.
'''

import pandas as pd
import argparse
import os

def make_uORFs_bed(args):
    uORFsString = ""
    for index, row in args.iterrows():
        uORFString = row.chromosome + "\t" + str(row.start) + "\t" + str(row.stop) + "\t" + row.uORFids + "\t0\t" + row.strand + "\n"
        uORFsString = uORFsString + uORFString
    return(uORFsString)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Converts uORF annotation.csv to bed format')
    parser.add_argument("--input_csv_filepath", help='Path to write \
                        merged csv output')
    parser.add_argument("--output_bed_filepath", help='Path to write \
                        merged bed6 output')
    args = parser.parse_args()
    uORFsdf = pd.read_csv(args.input_csv_filepath)
    uORFsbed = make_uORFs_bed(uORFsdf)
    f = open(args.output_bed_filepath, 'wt', encoding='utf-8')
    f.write(uORFsbed)


if __name__ == '__main__':
    main()
