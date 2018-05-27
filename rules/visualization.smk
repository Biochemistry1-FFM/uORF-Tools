rule genomeSamToolsIndex:
    input:
        rules.retrieveGenome.output
    output:
        "index/genomeSamtools/sizes.genome"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    params:
    log: "logs/_{condition}_{sampleid}.log"
    shell:
        "mkdir -p index/genomeSamtools, samtools faidx {rules.retrieveGenome.output}; cut -f1,2 genomes/genome.fa.fai > index/genomeSamtools/sizes.genome"

rule wig:
    input:
        expand("bam/{method}_{condition}_{sampleid}/Aligned.sortedByCoord.out.bam", method=config["methods"], condition=config["conditions"], sampleid=config["sampleids"]),
        rules.output.genomeSamToolsIndex
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
        rules.output.genomeSamToolsIndex
    output:
        "tracks/annotation.bb"
    conda:
        "../envs/bed.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; bedToBigBed {input[0]} {input[1]} tracks/annotation.bb"
