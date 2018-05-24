import os
import re
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

ADAPTERS=config["adapter"]
METHODS=config["methods"]
CONDITIONS=config["conditions"]
SAMPLEIDS=config["sampleids"]

rule all:
   input:
       expand("bam/{method}_{condition}_{sampleid}/Aligned.sortedByCoord.out.bam", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       expand("ribotaper/{condition}_{sampleid}/ORFs_max_filt", condition=CONDITIONS, sampleid=SAMPLEIDS),
       expand("uORFs/sfactors_lprot_{method}.csv", condition=CONDITIONS, method=METHODS),
       expand("uORFs/ncounts_lprot_{method}.csv", condition=CONDITIONS, method=METHODS)
onsuccess:
    print("Done, no error")

rule trim:
    input:
        expand("fastq/{method}_{condition}_{sampleid}.fastq.gz", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS)
    output:
        "trimmed/{method}_{condition}_{sampleid}.fastq"
    params: ada=ADAPTERS
    conda:
        "envs/cutadapt.yaml"
    threads: 20
    shell:
        "mkdir -p trimmed; cutadapt -a {params.ada} -j {threads} -u 1 -q 20 -O 1 -m 15 --trim-n -o {output} {input[0]}"

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
        "mkdir -p index/genomeStar; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir index/genomeStar --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100"

rule map:
    input:
        expand("norRNA/{method}_{condition}_{sampleid}.fastq", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
        rules.genomeIndex.output
    output:
        "bam/{method}_{condition}_{sampleid}/Aligned.sortedByCoord.out.bam"
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
        "mkdir -p ribotaper/ribotaper_annotation; create_annotations_files.bash {input[0]} {input[1]} true false ribotaper/ribotaper_annotation"

rule ribotaperMetaplot:
    input:
        rules.map.output,
        rules.ribotaperAnnotation.output
    output:
        "metaplots/{method}_{condition}_{sampleid}"
    conda:
        "envs/ribotaper.yaml"
    threads: 1
    shell:
        "mkdir -p ribotaper/metaplots; create_metaplots.bash {input[0]} {input[1]} {output[0]}"

rule ribotaper:
    input:
        fp=expand("bam/FP_{condition}_{sampleid}/Aligned.sortedByCoord.out.bam", condition=CONDITIONS, sampleid=SAMPLEIDS), total=expand("bam/Total_{condition}_{sampleid}/Aligned.sortedByCoord.out.bam", condition=CONDITIONS, sampleid=SAMPLEIDS),
        annotation=rules.ribotaperAnnotation.output
    output:
        "ribotaper/{condition}_{sampleid}/ORFs_max_filt",
        "ribotaper/{condition}_{sampleid}/Final_ORF_results.pdf"
    conda:
        "envs/ribotaper.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    log: "logs/ribotaper_{condition}_{sampleid}.log"
    shell:
        "mkdir -p {params.prefix}; cd {params.prefix}; Ribotaper.sh ../../{input.fp[0]} ../../{input.total[0]} ../../ribotaper/ribotaper_annotation/ 27,29,30,31 11,11,12,12 {threads}"

rule ribotaperMerge:
    input:
        ctrl=expand("ribotaper/ctrl_{sampleid}/ORFs_max_filt", sampleid=SAMPLEIDS), treat=expand("ribotaper/treat_{sampleid}/ORFs_max_filt", sampleid=SAMPLEIDS),
    output:
        "ribotaper/{sampleid}/Merged_uORF_results.csv"
    conda:
        "envs/uorftools.yaml"
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
       "bam/{method}_{condition}_{sampleid}/Aligned.sortedByCoord.out.bam"
    output:
        "bam/{method}_{condition}_{sampleid}.bam"
    threads: 1
    params:
        cwd=lambda wildcards, output: (os.getcwd())
    log: "logs/maplink.log"
    shell:
        "ln -s {params.cwd}/{input[0]} {params.cwd}/{output[0]};"

rule normalizedCounts:
    input:
        expand("bam/{method}_{condition}_{sampleid}.bam", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
        rules.longestTranscript.output
    output:
        "uORFs/ncounts_lprot_{method}.csv",
        "uORFs/sfactors_lprot_{method}.csv"
    conda:
        "envs/uorftools.yaml"
    threads: 20
    shell:
        "mkdir -p uORFs; uORF-Tools/generate_normalized_counts_longest_protein.R -r -b bam/ -a {input[1]} -s uORFs/sfactors_lprot_{method}.csv -n uORFs/ncounts_lprot_{method}.csv -t {method}"

# Import rules

include: "rules/rrnafiltering.smk"
