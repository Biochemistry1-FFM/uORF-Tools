import os
import re
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
ADAPTERS=config["ADAPTERS"]
RRNADB=config["RRNADBS"]
METHODS=["FP","Total"]
CONDITIONS=["ctrl","treat"]
SAMPLEIDS=["1-2"]


rule all:
   input:
       ["index/rRNA/{rrnadb}.bursttrie_0.dat".format(rrnadb=re.sub("-id\d+.fasta","",rrnadb)) for rrnadb in RRNADB],
       expand("trimmed/{method}_{condition}_{sampleid}.fastq", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
       expand("norRNA/{method}_{condition}_{sampleid}.fastq", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS)
       #"index/rRNA/{rrnadb}.kmer_0.dat".format(rrnadb=rrnadb) for rrnadb in RRNADB,
       #"index/rRNA/{rrnadb}.pos_0.dat".format(rrnadb=rrnadb) for rrnadb in RRNADB,
       #"index/rRNA/{rrnadb}.stats".format(rrnadb=rrnadb) for rrnadb in RRNADB]

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

rule rrnaretrieve:
    input:
        HTTP.remote("https://github.com/biocore/sortmerna/raw/master/rRNA_databases/{rrnadb}",keep_local=True,allow_redirects=True)
    output:
        "rRNA_databases/{rrnadb}"
    run:
        outputName = os.path.basename(input[0])
        shell("mkdir -p rRNA_databases; mv {input} rRNA_databases/{outputName}")

#define indexfiles function

def indexfiles (RRNADB):
    indexstring=""
    for rrnadb in RRNADB:
        rrnaprefix = rrnadb.replace(".fasta","")
        dbstring = "./rRNA_databases/" + rrnadb + ",./index/rRNA/" + rrnaprefix + "-db:"
        indexstring = dbstring + indexstring
    return str(indexstring)

rule rrnaindex:
    input:
        ["rRNA_databases/{rrnadb}".format(rrnadb=rrnadb) for rrnadb in RRNADB]
    output:
        ["index/rRNA/{rrnadb}.bursttrie_0.dat".format(rrnadb=re.sub("-id\d+.fasta","",rrnadb)) for rrnadb in RRNADB]
        #["index/rRNA/{rrnadb}.kmer_0.dat".format(rrnadb=rrnadb) for rrnadb in RRNADB]
        #["index/rRNA/{rrnadb}.pos_0.dat".format(rrnadb=rrnadb) for rrnadb in RRNADB]
        #["index/rRNA/{rrnadb}.stats".format(rrnadb=rrnadb) for rrnadb in RRNADB]
        #"index/rRNA/{rrnadb}.bursttrie_0.dat"
        #"index/rRNA/{rrnadb}.kmer_0.dat"
        #"index/rRNA/{rrnadb}.pos_0.dat"
        #"index/rRNA/{rrnadb}.stats"
    conda:
        "envs/sortmerna.yaml"
    shell:
        "mkdir -p index/rRNA; indexdb_rna --ref ./rRNA_databases/silva-euk-18s-id95.fasta,./index/rRNA/silva-euk-18s:./rRNA_databases/silva-euk-28s-id98.fasta,./index/rRNA/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rRNA/rfam-5s-database:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rRNA/rfam-5.8s-database"

rule rrnafilter:
    input:
        expand("trimmed/{method}_{condition}_{sampleid}.fastq", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS)
    output:
        "norRNA/{method}_{condition}_{sampleid}.fastq"
    conda:
        "envs/sortmerna.yaml"
    threads: 20
    shell:
        "mkdir -p norRNA; mkdir -p norRNA/rRNA; sortmerna -a {threads} --ref ./rRNA_databases/silva-euk-18s-id95.fasta,./index/rRNA/silva-euk-18s:./rRNA_databases/silva-euk-28s-id98.fasta,./index/rRNA/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rRNA/rfam-5s-database:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rRNA/rfam-5.8s-database --reads {input[0]} --num_alignments 1 --fastx --aligned norRNA/rRNA/reject --other {output}"
