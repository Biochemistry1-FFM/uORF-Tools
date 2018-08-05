#!/usr/bin/env python
'''This script creates a summary table containing information on the number of reads lost during each data prepration step of the workflow i.e. trimming, rRNA removal and mapping.'''


import argparse
import pandas as pd
import subprocess

# Create column names for data frame
def getName(row):
	name = "{}-{}-{}".format(row['method'], row['condition'], row['replicate'])
	return name

# Count number of raw reads
def countRaw(reads, folder):
	count = subprocess.Popen('zcat ' + folder + reads + ' | wc -l', shell = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = count.communicate()
	return stdout.decode('utf-8')

# Get file names of trimmed files
def getFileNames(folder):
	fileName = subprocess.Popen('for file in ' + folder + 'trimmed/*.fastq; do echo $file; done', shell = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = fileName.communicate()
	return stdout.decode('utf-8').split()

# Count number of trimmed reads, is also used for counting the number of no rRNA reads (both are uncompressed fastq files)
def countTrimmed(fileName):
	count = subprocess.Popen('wc -l ' + fileName, shell = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = count.communicate()
	return stdout.decode('utf-8').split()[0]

# Get file names of no rRNA reads
def getFileNamesNorRNA(folder):
	fileName = subprocess.Popen('for file in ' + folder + 'norRNA/*.fastq; do echo $file; done', shell = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = fileName.communicate()
	return stdout.decode('utf-8').split()

# Get file names of mapped reads
def getFileNamesMapped(folder):
	fileName = subprocess.Popen('for file in ' + folder + 'maplink/*.bam; do echo $file; done', shell = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = fileName.communicate()
	return stdout.decode('utf-8').split()

# Count number of mapped reads
def countMapped(fileName):
	count = subprocess.Popen('samtools view -c ' + fileName, shell = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = count.communicate()
	return stdout.decode('utf-8').split()[0]

# Create output table and calculate % losses
def createTable(sheet, folder):
	sheet = pd.read_table(sheet)
	names = sheet.apply(lambda row: getName(row), axis = 1)
	df = pd.DataFrame(columns = names, index = ['No. of raw reads', 'No. of trimmed reads', '% loss', 'No. of reads w/o rRNA', '% loss', 'No. of mapped reads', '% loss', '% total loss'])
	raw = [int(int(countRaw(i, folder).strip())/4) for i in sheet['fastqFile']]
	df.iloc[0] = raw 
	trimmedFiles = getFileNames(folder)
	trimmed = [int(int(countTrimmed(i))/4) for i in trimmedFiles]
	df.iloc[1] = trimmed
	norRNAFiles = getFileNamesNorRNA(folder)
	norRNA = [int(int(countTrimmed(i))/4) for i in norRNAFiles]
	df.iloc[3] = norRNA
	mappedFiles = getFileNamesMapped(folder)
	mapped = [int(countMapped(i)) for i in mappedFiles]
	df.iloc[5] = mapped
	df.iloc[2] = df.apply(lambda row: round(100 - (row.loc["No. of trimmed reads"] / row.loc["No. of raw reads"] * 100), 2))
	df.iloc[4] = df.apply(lambda row: round(100 - (row.loc["No. of reads w/o rRNA"] / row.loc["No. of trimmed reads"] * 100), 2))
	df.iloc[6] = df.apply(lambda row: round(100 - (row.loc["No. of mapped reads"] / row.loc["No. of reads w/o rRNA"] * 100), 2))
	df.iloc[7] = df.apply(lambda row: round(100 - (row.loc["No. of mapped reads"] / row.loc["No. of raw reads"] * 100), 2))
	return df

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Creates a summary table of read processing steps.')
    parser.add_argument('--sample_tsv', help = 'Path to sample sheet')
    parser.add_argument("--project_folder", help='Path to uORF-Tools project folder')
    parser.add_argument("--output_tsv_filepath", help='Path to write results table')
    args = parser.parse_args()
    df = createTable(args.sample_tsv, args.project_folder)
    df.to_csv(args.output_tsv_filepath, sep = "\t")


if __name__ == '__main__':
    main()