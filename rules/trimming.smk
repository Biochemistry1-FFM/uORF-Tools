def getfastq(wildcards):
    print(wildcards)
    return samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["fastqFile"]].dropna()

rule trim:
    input:
        reads=getfastq
    output:
        "trimmed/{method,[a-zA-Z]+}-{condition,[a-zA-Z]+}-{replicate,\d+}.fastq"
    params:
        ada=lambda wildcards, output: ("" if not ADAPTERS else (" -a " + ADAPTERS)),
        prefix=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.reads[0]))[0])[0])
    conda:
        "../envs/trimgalore.yaml"
    threads: 1
    shell:
        "mkdir -p trimmed; trim_galore {params.ada} --phred33 --output_dir trimmed/ --trim-n --suppress_warn --dont_gzip fastq/{params.prefix}.fastq.gz; mv trimmed/{params.prefix}_trimmed.fq {output}"

rule fastqcraw:
    input:
        "trimmed/{method}-{condition}-{replicate}.fastq"
    output:
        "fastqcraw/{method}-{condition}-{replicate}_raw.html"
    conda:
        "../envs/fastqc.yaml"
    threads: 6
    shell:
        "mkdir -p fastqcraw; fastqc -o raw -t {threads} {input}"

rule fastqctrimmed:
    input:
        "trimmed/{method}-{condition}-{replicate}.fastq"
    output:
        "fastqctrimmed/{method}-{condition}-{replicate}_trimmed.html"
    conda:
        "../envs/fastqc.yaml"
    threads: 6
    shell:
        "mkdir -p fastqctrimmed; fastqc -o trimmed -t {threads} {input}"
