rule genomeSize:
    input:
        rules.genomeSamToolsIndex.output
    output:
        "genomes/sizes.genome"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    params:
    log: "logs/genomeSamToolsIndex.log"
    shell:
        "mkdir -p genomes; cut -f1,2 {input[0]} > genomes/sizes.genome"


rule wig:
    input:
        expand("bam/{method}_{condition}_{sampleid}/Aligned.sortedByCoord.out.bam", method=config["methods"], condition=config["conditions"], sampleid=config["sampleids"]),
        rules.genomeSize.output
    output:
        "tracks/{method}_{condition}_{sampleid}.wig"
    conda:
        "../envs/wig.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    log: "logs/ribotaper_{condition}_{sampleid}.log"
    shell:
        "mkdir -p tracks; bam2wig.py -i {input[0]} -s {input[1]} -o tracks/{prefix}"

rule annotationBed:
    input:
        rules.retrieveAnnotation.output
    output:
        "tracks/annotation.bed"
    conda:
        "../envs/bed.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; gtf2bed < {input[0]} > tracks/annotation.bed"

rule annotationBigBed:
    input:
        rules.annotationBed.output,
        rules.genomeSize.output
    output:
        "tracks/annotation.bb"
    conda:
        "../envs/bed.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; bedToBigBed {input[0]} {input[1]} tracks/annotation.bb"
