rule genomeIndex:
    input:
        rules.retrieveGenome.output,
        rules.retrieveAnnotation.output
    output:
        "genomeStarIndex",
        #"genomeStarIndex/chrName.txt",
        #"genomeStarIndex/genomeParameters.txt"
    conda:
        "../envs/star.yaml"
    threads: 20
    params:
        indexpath=lambda wildcards: ("" if not INDEXPATH else (INDEXPATH))
    log:
        "logs/genomeIndex.log"
    shell:
        #"ln -T -s {params.indexpath} genomeStarIndex"
        #"if [ -d {params.indexpath} ]; then ln -T -s {params.indexpath} genomeStarIndex; echo \"Index linked\"; else echo \"Computing STAR index\"; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir genomeStarIndex --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100 2> {log}; fi"
        "if [ -d {params.indexpath} ]; then ln -T -s {params.indexpath} genomeStarIndex; echo \"Index linked\"; else mkdir -p genomeStarIndex; echo \"Computing STAR index\"; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir genomeStarIndex --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100 2> {log}; fi"

rule map:
    input:
        fastq="norRNA/{method}-{condition}-{replicate}.fastq",
        index=rules.genomeIndex.output
    output:
        "bam/{method, [a-zA-Z]+}-{condition, [a-zA-Z]+}-{replicate,\d+}/Aligned.sortedByCoord.out.bam"
    conda:
        "../envs/star.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    log:
        "logs/{method, [a-zA-Z]+}-{condition, [a-zA-Z]+}-{replicate,\d+}_star.log"
    shell:
        "mkdir -p bam; STAR --genomeDir genomeStarIndex --readFilesIn {input.fastq} --outFileNamePrefix {params.prefix}/ --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --runThreadN {threads} 2> {log}"

rule maplink:
    input:
        "bam/{method}-{condition}-{replicate}/Aligned.sortedByCoord.out.bam"
    output:
        "maplink/{method, [a-zA-Z]+}-{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink; ln -s {params.inlink} {params.outlink}"

rule ribomaplink:
    input:
        "bam/{method}-{condition}-{replicate}/Aligned.sortedByCoord.out.bam"
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
        "bam/{method}-{condition}-{replicate}/Aligned.sortedByCoord.out.bam"
    output:
        "maplink/{method, RNA}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNA/; ln -s {params.inlink} {params.outlink}"
