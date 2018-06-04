import os
import re
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

ADAPTERS=config["adapter"]
METHODS=config["methods"]
CONDITIONS=config["conditions"]
SAMPLEIDS=config["sampleids"]

(OUTWIGS) = glob_wildcards("fastq/{outwig}.fastq.gz")
print(OUTWIGS)

rule all:
   input:
       expand("fastqc/{method}-{condition}-{sampleid}-fastqc.html", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       expand("bam/{method}-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       expand("ribotaper/{condition}-{sampleid}/ORFs_max_filt", condition=CONDITIONS, sampleid=SAMPLEIDS),
       expand("bam/{method}-{condition}-{sampleid}.bam", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       expand("tracks/{outwig}.wig", outwig=OUTWIGS),
       "tracks/annotation.bb"
onsuccess:
    print("Done, no error")

rule trim:
    input:
        expand("fastq/{method}-{condition}-{sampleid}.fastq.gz", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS)
    output:
        "trimmed/{method}-{condition}-{sampleid}.fastq"
    params:
        ada=ADAPTERS,
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    conda:
        "envs/trimgalore.yaml"
    threads: 20
    prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    shell:
        "mkdir -p trimmed; trim_galore -a {params.ada} --phred33 --output_dir trimmed/ --trim-n --suppress_warn --dont_gzip {input[0]}; mv {params.prefix}_trimmed.fq {params.prefix}.fastq"

rule fastqc:
    input:
        rules.trim.output
    output:
        "fastqc/{method}-{condition}-{sampleid}_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    threads: 6
    shell:
        "mkdir -p fastqc; fastqc -o fastqc -t {threads} {input}"

rule retrieveGenome:
    input:
        "genome.fa"
    output:
        "genomes/genome.fa"
    threads: 20
    shell:
        "mkdir -p genomes; cp genome.fa genomes/"

rule retrieveAnnotation:
    input:
        "annotation.gtf"
    output:
        "annotation/annotation.gtf"
    threads: 20
    shell:
        "mkdir -p annotation; cp annotation.gtf annotation/"

rule genomeIndex:
    input:
        rules.retrieveGenome.output,
        rules.retrieveAnnotation.output
    output:
        "index/genomeStar/chrLength.txt",
        "index/genomeStar/chrName.txt",
        "index/genomeStar/genomeParameters.txt"
    conda:
        "envs/star.yaml"
    threads: 20
    shell:
        "mkdir -p index/genomeStar; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir index/genomeStar --genomeFastaFiles {input[0]}" #--sjdbGTFfile {input[1]} --sjdbOverhang 100"

ruleorder: map > maplink

rule map:
    input:
        expand("norRNA/{method}-{condition}-{sampleid}.fastq", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
        rules.genomeIndex.output
    output:
        "bam/{method}-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam"
    conda:
        "envs/star.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    shell:
        "mkdir -p bam; STAR --genomeDir index/genomeStar --readFilesIn {input[0]} --outFileNamePrefix {params.prefix}/ --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 --alignEndsType EndToEnd --runThreadN {threads}"

rule ribotaperAnnotation:
    input:
        rules.retrieveAnnotation.output,
        rules.retrieveGenome.output
    output:
        "ribotaper/ribotaper_annotation/start_stops_FAR.bed"
    conda:
        "envs/ribotaper.yaml"
    threads: 1
    shell:
        "mkdir -p ribotaper/ribotaper_annotation; create_annotations_files.bash {input[0]} {input[1]} false false ribotaper/ribotaper_annotation"

rule ribotaperMetaplot:
    input:
        rules.map.output,
        rules.ribotaperAnnotation.output
    output:
        "metaplots/{method}-{condition}-{sampleid}"
    conda:
        "envs/ribotaper.yaml"
    threads: 1
    shell:
        "mkdir -p ribotaper/metaplots; create_metaplots.bash {input[0]} {input[1]} {output[0]}"

rule genomeSamToolsIndex:
    input:
        rules.retrieveGenome.output
    output:
        "genomes/genome.fa.fai"
    conda:
        "envs/samtools.yaml"
    threads: 1
    params:
    shell:
        "samtools faidx {rules.retrieveGenome.output}"

rule ribotaper:
    input:
        fp=expand("bam/FP-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam", condition=CONDITIONS, sampleid=SAMPLEIDS), total=expand("bam/Total-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam", condition=CONDITIONS, sampleid=SAMPLEIDS),
        annotation=rules.ribotaperAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output
    output:
        "ribotaper/{condition}-{sampleid}/ORFs_max_filt",
        "ribotaper/{condition}-{sampleid}/Final_ORF_results.pdf"
    conda:
        "envs/ribotaper.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    shell:
        "mkdir -p {params.prefix}; cd {params.prefix}; Ribotaper.sh ../../{input.fp[0]} ../../{input.total[0]} ../../ribotaper/ribotaper_annotation/ 27,29,30,31 11,11,12,12 {threads}"

rule ribotaperMerge:
    input:
        ctrl=expand("ribotaper/ctrl-{sampleid}/ORFs_max_filt", sampleid=SAMPLEIDS), treat=expand("ribotaper/treat_{sampleid}/ORFs_max_filt", sampleid=SAMPLEIDS),
    output:
        "ribotaper/{sampleid}/Merged_uORF_results.csv"
    conda:
        "envs/uorftoolspython.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    shell:
        "mkdir -p {params.prefix}; uORF-Tools/scripts/ribotaper_merge_incl_length.py {input.ctrl} {input.treat} --output_csv_filepath {params.prefix}/Merged_uORF_results.csv"

rule longestTranscript:
    input:
        rules.retrieveAnnotation.output
    output:
        "uORFs/longest_protein_coding_transcripts.gtf"
    conda:
        "envs/uorftools.yaml"
    threads: 20
    shell:
        "mkdir -p uORFs; uORF-Tools/longest_orf_transcript.py -a {input} -o {output}"

rule maplink:
    input:
      expand("bam/{method}-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam", method=config["methods"], condition=config["conditions"], sampleid=config["sampleids"])
    output:
      expand("bam/{method}-{condition}-{sampleid}.bam", method=config["methods"], condition=config["conditions"], sampleid=config["sampleids"])
    params:
        cwd=os.getcwd()
    threads: 1
    run:
        for f in input:
                str=f
                outfile=str.replace("/Aligned.sortedByCoord.out.bam", ".bam")
                shell("ln -s {params.cwd}/{f} {params.cwd}/{outfile}")

rule normalizedCounts:
    input:
        rules.longestTranscript.output,
        rules.maplink.output
    output:
        "uORFs/sfactors_lprot_FP.csv",
        "uORFs/sfactors_lprot_Total.csv",
        "uORFs/ncounts_lprot_FP.csv",
        "uORFs/ncounts_lprot_Total.csv"
    conda:
        "envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/generate_normalized_counts_longest_protein.R -r -b bam/ -a {input[0]} -s uORFs/sfactors_lprot_FP.csv -n uORFs/ncounts_lprot_FP.csv -t FP;  uORF-Tools/generate_normalized_counts_longest_protein.R -r -b bam/ -a {input[0]} -s uORFs/sfactors_lprot_Total.csv -n uORFs/ncounts_lprot_Total.csv -t Total")

# Import rules

include: "rules/rrnafiltering.smk"
include: "rules/visualization.smk"
