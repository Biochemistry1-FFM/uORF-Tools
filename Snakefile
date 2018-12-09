import os
import re
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.3.1")

INDEXPATH=config["genomeindexpath"]
UORFANNOTATIONPATH=config["uorfannotationpath"]
CODONS=config["alternativestartcodons"]

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index(["method", "condition", "replicate"], drop=False)
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(samples)
#samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")
report: "report/workflow.rst"

rule all:
   input:
       expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples()),
       expand("ribotish/{sample.condition}-{sample.replicate}-newORFs.tsv_all.txt", sample=samples.itertuples()),
       expand("tracks/{sample.method}-{sample.condition}-{sample.replicate}.bw", sample=samples.itertuples()),
       "uORFs/merged_uORFs.csv",
       "tracks/annotation.bb",
       "uORFs/sfactors_lprot.csv",
       "uORFs/xtail_uORFs.csv",
       "uORFs/xtail_cds.csv",
       "uORFs/uORFs_regulation.tsv",

onsuccess:
    print("Done, no error")

#preprocessing
include: "rules/preprocessing.smk"
#bootstrap
include: "rules/bootstrap.smk"
#Visualization
include: "rules/visualization.smk"
#Ribotish
include: "rules/ribotish.smk"
#uORF-tools
include: "rules/uORF-Tools.smk"
