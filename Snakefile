import os
import re
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.3.1")

INDEXPATH = config["genomeindexpath"]
UORFANNOTATIONPATH = config["uorfannotationpath"]
CODONS = config["alternativestartcodons"]

samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])


def replicate_check(samples):
    if(samples['replicate'].nunique() < 2):
        print("Warning: Please make sure your experiment contains replicates!")


with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(samples)

validate(samples, schema="schemas/samples.schema.yaml")

onstart:
    if not os.path.exists("logs"):
        os.makedirs("logs")
    replicate_check(samples)

report: "report/workflow.rst"

rule all:
    input:
        expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples()),
        expand("ribotish/{sample.condition}-{sample.replicate}-newORFs.tsv_all.txt", sample=samples.itertuples()),
        expand("tracks/{sample.method}-{sample.condition}-{sample.replicate}.bw", sample=samples.itertuples()),
        "uORFs/merged_uORFs.csv",
        "tracks/annotation.bb",
        "uORFs/sfactors_lprot.csv",
        "uORFs/uORFs_regulation.tsv",

onsuccess:
    print("Done, no error")

# Preprocessing
include: "rules/preprocessing.smk"
# Bootstrap
include: "rules/bootstrap.smk"
# Visualization
include: "rules/visualization.smk"
# Ribotish
include: "rules/ribotish.smk"
# uORF-tools
include: "rules/uORF-Tools.smk"
