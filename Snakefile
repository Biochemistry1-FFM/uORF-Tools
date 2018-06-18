import os
import re
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

ADAPTERS=config["adapter"]

onstart:
   if not os.path.exists("logs"):
    os.makedirs("logs")

samples = pd.read_table(config["samples"], dtype=str).set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index
validate(samples, schema="schemas/samples.schema.yaml")

rule all:
   input:
       expand("trimmed/{method}-{condition}-{replicate}.fastq", **samples)
       #expand("fastqc/{method}-{condition}-{sampleid}_fastqc.html", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       #expand("bam/{method}-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       #expand("ribotaper/{condition}-{sampleid}/ORFs_max_filt", condition=CONDITIONS, sampleid=SAMPLEIDS),
       #expand("bam/{method}-{condition}-{sampleid}.bam", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       #expand("tracks/{method}-{condition}-{sampleid}.wig", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       #"tracks/annotation.bb",
       #expand("metaplots/{method}-{condition}-{sampleid}", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS)
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
#include: "rules/rrnafiltering.smk"
#include: "rules/mapping.smk"

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

#include: "rules/ribotaper.smk"
#include: "rules/uORF-Tools.smk"
#include: "rules/visualization.smk"

