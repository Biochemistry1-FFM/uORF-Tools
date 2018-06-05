rule trim:
    input:
        expand("fastq/{method}-{condition}-{sampleid}.fastq.gz", method=METHODS, condition=CONDITIONS, sampleid=SAMPLEIDS)
    output:
        "trimmed/{method}-{condition}-{sampleid}.fastq"
    params:
        ada=lambda wildcards, output: ("" if not ADAPTERS else (" -a " + ADAPTERS)),
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    conda:
        "../envs/trimgalore.yaml"
    threads: 1
    shell:
        "mkdir -p trimmed; trim_galore {params.ada} --phred33 --output_dir trimmed/ --trim-n --suppress_warn --dont_gzip fastq/{params.prefix}.fastq.gz; mv trimmed/{params.prefix}_trimmed.fq trimmed/{params.prefix}.fastq"
