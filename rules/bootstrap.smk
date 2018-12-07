def getbam(wildcards):
    return tamples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["inputFile"]].dropna()

rule rawbamlink:
    input:
        bams=getbam
    output:
        "maplink/{method}-{condition}-{replicate}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + (os.path.splitext(input.bams[0])[0]) + ".bam"),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink; ln -s {params.inlink} {params.outlink}"

rule ribomaplink:
    input:
        "maplink/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method, RIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RIBO/; ln -s {params.inlink} {params.outlink}"

rule rnamaplink:
    input:
        "maplink/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method, RNA}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNA/; ln -s {params.inlink} {params.outlink}"
