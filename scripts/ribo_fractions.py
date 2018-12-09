#!/usr/bin/env python
'''Computes differences between riboseq reads between conditions
'''

import pandas as pd
import re
import argparse
import numpy as np
import os

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Computes differences between riboseq reads between conditions')
    parser.add_argument("--uORF_reads", help='Path to input file with uORF reads')
    parser.add_argument("--ORF_reads", help='Path to input file with ORF reads')
    args = parser.parse_args()

if __name__ == '__main__':
    main()
