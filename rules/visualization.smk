rule genomeSamToolsIndex:
    input:
        rules.retrieveGenome.output
    output:
        "index/genomeSamtools/sizes.genome"
    conda:
        "envs/samtools.yaml"
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
        "wig/{method}_{condition}_{sampleid}.wig"
    conda:
        "envs/wig.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    log: "logs/ribotaper_{condition}_{sampleid}.log"
    shell:
        "mkdir -p wig; bam2wig.py -i {input[0]} -s {input[1]} -o prefix"
