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
        "bam/{outwig}/Aligned.sortedByCoord.out.bam",
        rules.genomeSize.output
    output:
        "tracks/{outwig}.wig"
    conda:
        "../envs/wig.yaml"
    threads: 1
    log: "logs/wig.log"
    shell:
        "mkdir -p tracks; bam2wig.py -i {input[0]} -s {input[1]} -o tracks/{outwig}"

rule annotationBed:
    input:
        rules.retrieveAnnotation.output
    output:
        "tracks/annotation.bed"
    conda:
        "../envs/bed.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input[0]} | grep -v '\tgene\t' > tracks/annotation-woGenes.gtf; gtf2bed < tracks/annotation-woGenes.gtf > tracks/annotation.bed"

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
        "mkdir -p tracks; bedToBigBed -tab {input[0]} {input[1]} tracks/annotation.bb"
