def getbam(wildcards):
    return samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["inputFile"]].dropna()

rule rawbamlink:
    input:
        bams=getbam
    output:
        "maplink/{method}-{condition}-{replicate}.bam"
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.bams[0]))[0])[0])
    threads: 1
    shell:
        "mkdir -p maplink; cp bam/{params.prefix}.bam {output}"

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
