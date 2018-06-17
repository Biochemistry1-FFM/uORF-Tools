import os
import re
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

ADAPTERS=config["adapter"]
#METHODS=config["methods"]
#CONDITIONS=config["conditions"]
#SAMPLEIDS=config["sampleids"]

onstart:
   if not os.path.exists("logs"):
    os.makedirs("logs")

methods = pd.read_table(config["methods"]).set_index("method", drop=False)

conditions = pd.read_table(config["conditions"]).set_index("condition", drop=False)

samples = pd.read_table(config["samples"], dtype=str).set_index(["methods", "sample", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index


rule all:
   input:
       expand("fastqc/{method}-{condition}-{sampleid}_fastqc.html", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       #expand("bam/{method}-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       #expand("ribotaper/{condition}-{sampleid}/ORFs_max_filt", condition=CONDITIONS, sampleid=SAMPLEIDS),
       #expand("bam/{method}-{condition}-{sampleid}.bam", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       #expand("tracks/{method}-{condition}-{sampleid}.wig", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       #"tracks/annotation.bb"
onsuccess:
    print("Done, no error")

rule retrieveGenome:
    input:
        "genome.fa"
    output:
        "genomes/genome.fa"
    threads: 1
    shell:
        "mkdir -p genomes; cp genome.fa genomes/"

rule retrieveAnnotation:
    input:
        "annotation.gtf"
    output:
        "annotation/annotation.gtf"
    threads: 1
    shell:
        "mkdir -p annotation; cp annotation.gtf annotation/"

include: "rules/trimming.smk"
include: "rules/rrnafiltering.smk"
include: "rules/mapping.smk"

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

include: "rules/ribotaper.smk"

rule ribotaperMerge:
    input:
        ctrl=expand("ribotaper/ctrl-{sampleid}/ORFs_max_filt", sampleid=SAMPLEIDS), treat=expand("ribotaper/treat_{sampleid}/ORFs_max_filt", sampleid=SAMPLEIDS),
    output:
        "ribotaper/{sampleid}/Merged_uORF_results.csv"
    conda:
        "envs/uorftoolspython.yaml"
    threads: 6
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
        "mkdir -p uORFs; uORF-Tools/scripts/longest_orf_transcript.py -a {input} -o {output}"

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
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_normalized_counts_longest_protein.R -r -b bam/ -a {input[0]} -s uORFs/sfactors_lprot_FP.csv -n uORFs/ncounts_lprot_FP.csv -t FP;  uORF-Tools/scripts/generate_normalized_counts_longest_protein.R -r -b bam/ -a {input[0]} -s uORFs/sfactors_lprot_Total.csv -n uORFs/ncounts_lprot_Total.csv -t Total")

include: "rules/visualization.smk"
