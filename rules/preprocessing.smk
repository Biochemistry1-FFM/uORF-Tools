rule retrieveGenome:
    input:
        "genome.fa"
    output:
        "genomes/genome.fa"
    threads: 1
    shell:
        "mkdir -p genomes; cp genome.fa genomes/"

rule retrieveAnnotation:
    input:
        "annotation.gtf"
    output:
        "annotation/annotation.gtf"
    threads: 1
    shell:
        "mkdir -p annotation; cp annotation.gtf annotation/"

