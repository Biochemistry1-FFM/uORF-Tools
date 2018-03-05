#!/usr/bin/env python
import re
import argparse


def gtf_gene_id(gtf):
    r"""Get gene id from GTF line
    >>> gtf_gene_id('chr1\tENSEMBL\tUTR\t70006\t70108\t.\t+\t.\tgene_id \"ENSG00000186092.5\"; transcript_id "ENST00000335137.4"; gene_type "protein_coding"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_name "OR4F5-201"; exon_number 1; exon_id "ENSE00002319515.2"; level 3; protein_id "ENSP00000334393.3"; transcript_support_level "NA"; tag "basic"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.3";')
    'ENSG00000186092.5'
    """
    gtffields = gtf.split("\t")
    extravalues = gtffields[8].split(";")
    gene_id = extravalues[0]
    line_gene = re.findall(r'gene_id', gene_id)
    if line_gene:
        dkey = re.sub('gene_id "', "", gene_id)
        dkey2 = re.sub('"', "", dkey)
        return dkey2
    else:
        return None


def gtf_transkript_length(orf):
    r"""Get transkript length from GTF
    >>> gtf_transkript_length("chr1\tENSEMBL\tUTR\t70006\t70108\t.\t+\t.\tgene_id \"ENSG00000186092.5\"; transcript_id \"ENST00000335137.4\"; gene_type \"protein_coding\"; gene_name \"OR4F5\"; transcript_type \"protein_coding\"; transcript_name \"OR4F5-201\"; exon_number 1; exon_id \"ENSE00002319515.2\"; level 3; protein_id \"ENSP00000334393.3\"; transcript_support_level \"NA\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS30547.1\"; havana_gene \"OTTHUMG00000001094.3\";")
    102
    >>> gtf_transkript_length("chr1\tENSEMBL\tUTR\t70108\t70006\t.\t-\t.\tgene_id \"ENSG00000186092.5\"; transcript_id \"ENST00000335137.4\"; gene_type \"protein_coding\"; gene_name \"OR4F5\"; transcript_type \"protein_coding\"; transcript_name \"OR4F5-201\"; exon_number 1; exon_id \"ENSE00002319515.2\"; level 3; protein_id \"ENSP00000334393.3\"; transcript_support_level \"NA\"; tag \"basic\"; tag \"appris_principal_1\"; tag \"CCDS\"; ccdsid \"CCDS30547.1\"; havana_gene \"OTTHUMG00000001094.3\";")
    102
    """
    gtffields = orf.split("\t")
    strand = gtffields[6]
    [strand] = re.findall(r'^[-+]$', strand)
    start = gtffields[3]
    end = gtffields[4]
    if strand == "+":
        length = int(end) - int(start)
        return length
    else:
        rlength = int(start) - int(end)
        return rlength


def main():
    parser = argparse.ArgumentParser(description='Extract longest transcript per gene id from annotation gtf.')
    parser.add_argument("-a","--annotation-gtf-filepath", help='Path to annotation gtf file')
    parser.add_argument("-o","--output-gtf-filepath", help='Path to output gtf file')
    args = parser.parse_args()
    orfs = []
    with open(args.annotation-gtf-filepath) as origin_file:
        for line in origin_file:
            lineCoding = re.findall(r'protein_coding', line)
            if lineCoding:
                gtffields = line.split("\t")
                if gtffields[2] == "transcript":
                    orfs.append(line)
    transcripts = {}
    for orf in orfs:
        dkey = gtf_gene_id(orf)
        if dkey is not None:
            if dkey in transcripts:
                transcripts[dkey].extend([orf])
            else:
                transcripts[dkey] = [orf]
                
    outfile = open("args.output-gtf-filepath","w+")
    for key, value in transcripts.items():
        max_orf = ""
        max_length = 0
        for orf in value:
            orf_length = gtf_transkript_length(orf)
            if orf_length > max_length:
                max_orf = orf
                max_length = orf_length
        if max_orf:
            outfile.write(max_orf)


if __name__ == '__main__':
    main()
