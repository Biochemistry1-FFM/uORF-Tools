import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
RRNADB=config["RRNADBS"]

rule all:
   input:
       ["rRNA_database/{rrnadb}".format(rrnadb=rrnadb) for rrnadb in RRNADB]

rule rrnaretrieve:
    input:
        HTTP.remote("https://github.com/biocore/sortmerna/raw/master/rRNA_databases/{rrnadb}",keep_local=True,allow_redirects=True)
    output:
        "rRNA_database/{rrnadb}"
    run:
        outputName = os.path.basename(input[0])
        shell("mkdir -p rRNA_database; mv {input} rRNA_database/{outputName}")

#define indexfiles function

def indexfiles (RRNADB):
    indexstring=""
    for rrnadb in RRNADB:
        rrnaprefix = rrnadb.replace(".fasta","")
        dbstring = "./rRNA_databases/" + rrnadb + ",./index/rRNA/" + rrnaprefix + "-db:"
        indexstring = dbstring + indexstring
    return indexstring

rule rrnaindex:
    input:
        ["rRNA_database/{rrnadb}".format(rrnadb=rrnadb) for rrnadb in RRNADB]
    output:
        "index/rRNA/{rrnadb}.bursttrie_0.dat"
        "index/rRNA/{rrnadb}.kmer_0.dat"
        "index/rRNA/{rrnadb}.pos_0.dat"
        "index/rRNA/{rrnadb}.stats"
    conda:
        "envs/sortmerna.yaml"
    shell:
        "mkdir -p index/rRNA; indexdb_rna -a 20 --ref {indexfiles}"
