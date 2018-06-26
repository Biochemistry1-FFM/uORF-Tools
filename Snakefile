import os
import re
import pandas as pd
from snakemake.utils import validate, min_version
min_version("5.1.2")

ADAPTERS=config["adapter"]
INDEXPATH=config["genomeindexpath"]

onstart:
   if not os.path.exists("logs"):
    os.makedirs("logs")

samples = pd.read_table(config["samples"], dtype=str).set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index
validate(samples, schema="schemas/samples.schema.yaml")

#wildcard_constraints:
#   method="\[a-zA-Z]+",
#   condition="\[a-zA-Z]+",
#   replicate="\d+"


rule all:
   input:
       expand("fastqcraw/{method}-{condition}-{replicate}_raw.html", **samples),
       expand("fastqctrimmed/{method}-{condition}-{replicate}_trimmed.html", **samples),
       expand("fastqcrrnafilter/{method}-{condition}-{replicate}_rrnafilter.html", **samples),
       expand("ribotaper/{condition}-{replicate}/ORFs_max_filt", **samples),
       expand("tracks/{method}-{condition}-{replicate}.wig", **samples),
       "tracks/annotation.bb"
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
#ribotaper
include: "rules/ribotaper.smk"
#uORF-tools
include: "rules/uORF-Tools.smk"
#Visualization
include: "rules/visualization.smk"
