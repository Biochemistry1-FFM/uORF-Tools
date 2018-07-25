rule genomeSamToolsIndex:
    input:
        rules.retrieveGenome.output
    output:
        "genomes/genome.fa.fai"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    shell:
        "samtools faidx {rules.retrieveGenome.output}"

rule ribotish:
    input:
        fp="maplink/RIBO/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation=rules.retrieveAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
	bamindex=rules.genomeSize.output
    output:
        report=report("ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-newORFs.tsv", caption="../report/ribotaper.rst", category="Ribotish")
    conda:
        "../envs/ribotish.yaml"
    threads: 1
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotish.log"
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output.raw))
    shell:
        "mkdir -p ribotish; ribotish predict --longest -b {input.fp} -g {input.annotation} -f {input.genome} -o {output.report} 2> {log}"
