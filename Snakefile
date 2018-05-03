import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

#def _rrnadbs():
    #ret = list()
    #for x in list(config["RRNADBS"]):
    #    ret.append(str(x))

RRNADB=config["RRNADBS"]
#print(RRNADB[0])
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
