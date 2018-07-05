rule ribotaperMerge:
    input:
        expand("ribotaper/{sample.condition}-{sample.replicate}-newORFs.tsv", sample=samples.itertuples())
    output:
        "ribotaper/Merged_uORF_results.csv"
    conda:
        "../envs/uorftoolspython.yaml"
    threads: 6
    shell:
        "uORF-Tools/scripts/ribotaper_merge_incl_length.py {input} --max_length 400 --output_csv_filepath ribotaper/Merged_uORF_results.csv"

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
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_normalized_counts_longest_protein.R -r -b maplink -a {input[0]} -s uORFs/sfactors_lprot.csv;")

