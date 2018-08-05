import os
import re
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.2.1")

ADAPTERS=config["adapter"]
INDEXPATH=config["genomeindexpath"]
UORFANNOTATIONPATH=config["uorfannotationpath"]


onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_table(config["samples"], dtype=str).set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")
report: "report/workflow.rst"

rule all:
   input:
       expand("fastqc/raw/{sample.method}-{sample.condition}-{sample.replicate}-raw.html", sample=samples.itertuples()),
       expand("fastqc/trimmed/{sample.method}-{sample.condition}-{sample.replicate}-trimmed.html", sample=samples.itertuples()),
       expand("fastqc/norRNA/{sample.method}-{sample.condition}-{sample.replicate}-norRNA.html", sample=samples.itertuples()),
       expand("ribotish/{sample.condition}-{sample.replicate}-newORFs.tsv", sample=samples.itertuples()),
       expand("report/{sample.condition}-{sample.replicate}-qual.jpg", sample=samples.itertuples()),
       expand("tracks/{sample.method}-{sample.condition}-{sample.replicate}.bw", sample=samples.itertuples()),
       "uORFs/Merged_uORF_results.csv",
       "tracks/annotation.bb",
       "uORFs/sfactors_lprot.csv",
       "uORFs/xtail_uORFs.csv",
       "uORFs/xtail_cds.csv",
       "uORFs/uORF_regulation.tsv",
       "report/xtail_cds_fc.jpg",
       "report/xtail_cds_r.jpg",
       "report/xtail_cds_fc.jpg",
       "report/xtail_cds_r.jpg",
       "uORFs/summary_results.tsv"

onsuccess:
    print("Done, no error")

#Copying of user provided genome annotion
include: "rules/preprocessing.smk"
#Adaper removal and quality control
include: "rules/trimming.smk"
#removal of reads mapping to ribosomal rna genes
include: "rules/rrnafiltering.smk"
#mapping
include: "rules/mapping.smk"
#Visualization
include: "rules/visualization.smk"
#Ribotish
include: "rules/ribotish.smk"
#uORF-tools
include: "rules/uORF-Tools.smk"
#Report
include: "rules/report.smk"
