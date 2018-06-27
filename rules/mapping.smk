rule genomeIndex:
    input:
        rules.retrieveGenome.output,
        rules.retrieveAnnotation.output
    output:
        "index/genomeStar/chrLength.txt",
        "index/genomeStar/chrName.txt",
        "index/genomeStar/genomeParameters.txt"
    conda:
        "../envs/star.yaml"
    threads: 20
    params:
        indexpath=lambda wildcards: ("" if not INDEXPATH else (INDEXPATH))
    shell:
        "if [ -d {indexpath} ]; then ln -s {indexpath} index/genomeStar; else mkdir -p index/genomeStar; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir index/genomeStar --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100; fi"

#ruleorder: map > maplink

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
    shell:
        "mkdir -p bam; STAR --genomeDir index/genomeStar --readFilesIn {input.fastq} --outFileNamePrefix {params.prefix}/ --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --runThreadN {threads}"

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
