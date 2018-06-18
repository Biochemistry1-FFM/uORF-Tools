rule ribotaperAnnotation:
    input:
        annotation=rules.retrieveAnnotation.output,
        genome=rules.retrieveGenome.output
    output:
        "ribotaper/ribotaper_annotation/start_stops_FAR.bed"
    conda:
        "../envs/ribotaper.yaml"
    threads: 1
    shell:
        "mkdir -p ribotaper/ribotaper_annotation; create_annotations_files.bash {input.annotation} {input.genome} true false ribotaper/ribotaper_annotation"

rule ribotaperMetaplot:
    input:
        map=rules.map.output,
        annotation=rules.ribotaperAnnotation.output
    output:
        "metaplots/{method}-{condition}-{replicate}"
    conda:
        "../envs/ribotaper.yaml"
    threads: 1
    shell:
        "mkdir -p ribotaper/metaplots; create_metaplots.bash {input.map} {input.annotation} {output[0]}"

rule genomeSamToolsIndex:
    input:
        rules.retrieveGenome.output
    output:
        "genomes/genome.fa.fai"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    params:
    shell:
        "samtools faidx {rules.retrieveGenome.output}"

rule ribotaper:
    input:
        fp=expand("bam/RIBO-{condition}-{replicate}/Aligned.sortedByCoord.out.bam", **samples), total=expand("bam/RNA-{condition}-{replicate}/Aligned.sortedByCoord.out.bam", **samples),
        annotation=rules.ribotaperAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output
    output:
        "ribotaper/{condition}-{replicate}/ORFs_max_filt",
        "ribotaper/{condition}-{replicate}/Final_ORF_results.pdf"
    conda:
        "../envs/ribotaper.yaml"
    threads: 6
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    shell:
        "mkdir -p {params.prefix}; cd {params.prefix}; Ribotaper.sh ../../{input.fp[0]} ../../{input.total[0]} ../../ribotaper/ribotaper_annotation/ 27,29,30,31 11,11,12,12 {threads}"
