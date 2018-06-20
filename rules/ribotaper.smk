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
        map="maplink/{method}-{condition}-{replicate}.bam",
        annotation=rules.ribotaperAnnotation.output
    output:
        "metaplots/{method}-{condition}-{replicate}.plot"
    conda:
        "../envs/ribotaper.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p metaplots; mkdir -p ribotaper/metaplots/{params.prefix}; cd ribotaper/metaplots/{params.prefix}; create_metaplots.bash ../../../{input.map} ../../../{input.annotation} {params.prefix}; mv {params.prefix} ../../../{output}"

rule genomeSamToolsIndex:
    input:
        rules.retrieveGenome.output
    output:
        "genomes/genome.fa.fai"
    conda:
        "../envs/samtools.yaml"
    threads: 20
    params:
    shell:
        "samtools faidx -@ {threads} {rules.retrieveGenome.output}"

rule psiteOffset:
    input:
        mplot="metaplots/{method}-{condition}-{replicate}.plot"
    output:
        "offsets/{method, RIBO}-{condition}-{replicate}.offset"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell:
        "mkdir -p offsets; uORF-Tools/scripts/calculate_p_site_offset.R -i {input.mplot} -o {output}"


rule ribotaper:
    input:
        fp=expand("maplink/RIBO-{condition}-{replicate}.bam", **samples), 
        total=expand("maplink/RNA-{condition}-{replicate}.bam", **samples),
        offset=expand("offsets/RIBO-{condition}-{replicate}.offset", **samples),
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
        "mkdir -p {params.prefix}; offset={{$}}(<{input.offset}) cd {params.prefix}; Ribotaper.sh ../../{input.fp[0]} ../../{input.total[0]} ../../ribotaper/ribotaper_annotation/ {{$offset}} {threads}"
