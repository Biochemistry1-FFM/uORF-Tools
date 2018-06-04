rule genomeIndex:
    input:
        rules.retrieveGenome.output,
        rules.retrieveAnnotation.output
    output:
        "index/genomeStar/chrLength.txt",
        "index/genomeStar/chrName.txt",
        "index/genomeStar/genomeParameters.txt"
    conda:
        "envs/star.yaml"
    threads: 20
    shell:
        "mkdir -p index/genomeStar; STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir index/genomeStar --genomeFastaFiles {input[0]}" #--sjdbGTFfile {input[1]} --sjdbOverhang 100"

ruleorder: map > maplink

rule map:
    input:
        expand("norRNA/{method}-{condition}-{sampleid}.fastq", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS),
        rules.genomeIndex.output
    output:
        "bam/{method}-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam"
    conda:
        "envs/star.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    shell:
        "mkdir -p bam; STAR --genomeDir index/genomeStar --readFilesIn {input[0]} --outFileNamePrefix {params.prefix}/ --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 --alignEndsType EndToEnd --runThreadN {threads}"

