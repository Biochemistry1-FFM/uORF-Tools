def getbam(wildcards):
    return samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["inputFile"]].dropna()

rule rawbamlink:
    input:
        bam=getbam
    output:
        "maplink/{method}-{condition}-{replicate}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    conda:
        "../envs/trimgalore.yaml"
    threads: 1
    shell:
        "mkdir -p maplink; ln -s {params.inlink} {params.outlink}"
