rule ribotaperAnnotation:
    input:
        rules.retrieveAnnotation.output,
        rules.retrieveGenome.output
    output:
        "ribotaper/ribotaper_annotation/start_stops_FAR.bed"
    conda:
        "../envs/ribotaper.yaml"
    threads: 1
    shell:
        "mkdir -p ribotaper/ribotaper_annotation; create_annotations_files.bash {input[0]} {input[1]} true false ribotaper/ribotaper_annotation"

rule ribotaperMetaplot:
    input:
        rules.map.output,
        rules.ribotaperAnnotation.output
    output:
        "metaplots/{method}-{condition}-{sampleid}"
    conda:
        "../envs/ribotaper.yaml"
    threads: 1
    shell:
        "mkdir -p ribotaper/metaplots; create_metaplots.bash {input[0]} {input[1]} {output[0]}"

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
        fp=expand("bam/FP-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam", condition=CONDITIONS, sampleid=SAMPLEIDS), total=expand("bam/Total-{condition}-{sampleid}/Aligned.sortedByCoord.out.bam", condition=CONDITIONS, sampleid=SAMPLEIDS),
        annotation=rules.ribotaperAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output
    output:
        "ribotaper/{condition}-{sampleid}/ORFs_max_filt",
        "ribotaper/{condition}-{sampleid}/Final_ORF_results.pdf"
    conda:
        "../envs/ribotaper.yaml"
    threads: 6
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    shell:
        "mkdir -p {params.prefix}; cd {params.prefix}; Ribotaper.sh ../../{input.fp[0]} ../../{input.total[0]} ../../ribotaper/ribotaper_annotation/ 27,29,30,31 11,11,12,12 {threads}"
