rule ribobamindexlink:
    input:
        "maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        "maplink/{method, RIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.bai"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    log: "logs/{method}-{condition}-{replicate}_ribobamindexlink.log"
    shell:
        "mkdir -p maplink/RIBO/; ln -s {params.inlink} {params.outlink} 2> {log}"

rule ribotishQuality:
    input:
        fp="maplink/RIBO/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation="uORFs/longest_protein_coding_transcripts.gtf",
        samindex=rules.genomeSamToolsIndex.output,
        bamindex="maplink/RIBO/{condition}-{replicate}.bam.bai"
	#bamindex=expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam.bai", sample=samples.itertuples())
    output:
        reportpdf="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.pdf",
        reporttxt=report("ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.txt", caption="../report/ribotishquality.rst", category="Ribotish"),
	offsetparameters="maplink/RIBO/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.para.py"
    conda:
        "../envs/ribotish.yaml"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotishquality.log"
    shell:
        "mkdir -p ribotish; ribotish quality -p {threads} -b {input.fp} -g {input.annotation} -o {output.reporttxt} -f {output.reportpdf} 2> {log}"

rule ribotish:
    input:
        fp="maplink/RIBO/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation="uORFs/longest_protein_coding_transcripts.gtf",
        samindex=rules.genomeSamToolsIndex.output,
        bamindex="maplink/RIBO/{condition}-{replicate}.bam.bai",
        offsetparameters="maplink/RIBO/{condition}-{replicate}.bam.para.py"
    output:
        report=report("ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-newORFs.tsv_all.txt", caption="../report/ribotish.rst", category="Ribotish"),
        filtered="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-newORFs.tsv"
    params:
        report="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt",
        codons= lambda wildcards: ("" if not CODONS else (" --alt --altcodons " + CODONS))
    conda:
        "../envs/ribotish.yaml"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotish.log"
    shell:
        "mkdir -p ribotish; if grep -q \"offdict = {{'m0': {{}}}}\" {input.offsetparameters}; then mv {input.offsetparameters} {input.offsetparameters}.unused; fi; ribotish predict --longest -v {params.codons} -p {threads} -b {input.fp} -g {input.annotation} -f {input.genome} -o {output.filtered} 2> {log}"
