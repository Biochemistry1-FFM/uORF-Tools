rule riboMerge:
    input:
        expand("ribotish/{sample.condition}-{sample.replicate}-newORFs.tsv", sample=samples.itertuples())
    output:
        "uORFs/Merged_uORF_results.bed",
        "uORFs/Merged_uORF_results.csv"
    conda:
        "../envs/uorftoolspython.yaml"
    threads: 1
    shell:
        "mkdir -p uORFs; uORF-Tools/scripts/ribo_merge.py {input} --min_length 1 --max_length 400 --output_csv_filepath uORFs/Merged_uORF_results.csv --output_bed_filepath uORFs/Merged_uORF_results.bed"

rule longestTranscript:
    input:
        rules.retrieveAnnotation.output
    output:
        "uORFs/longest_protein_coding_transcripts.gtf"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell:
        "mkdir -p uORFs; uORF-Tools/scripts/longest_orf_transcript.py -a {input} -o {output}"

rule sizeFactors:
    input:
        rules.longestTranscript.output,
        expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples())
    output:
        "uORFs/sfactors_lprot.csv"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_size_factors.R -t uORF-Tools/samples.tsv -b maplink/ -a {input[0]} -s uORFs/sfactors_lprot.csv;")

