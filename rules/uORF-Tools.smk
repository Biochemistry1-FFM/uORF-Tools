rule ribotaperMerge:
    input:
        ctrl=expand("ribotaper/ctrl-{replicate}/ORFs_max_filt", **samples), treat=expand("ribotaper/treat_{replicate}/ORFs_max_filt", **samples),
    output:
        "ribotaper/{replicate}/Merged_uORF_results.csv"
    conda:
        "envs/uorftoolspython.yaml"
    threads: 6
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    shell:
        "mkdir -p {params.prefix}; uORF-Tools/scripts/ribotaper_merge_incl_length.py {input.ctrl} {input.treat} --output_csv_filepath {params.prefix}/Merged_uORF_results.csv"

rule longestTranscript:
    input:
        rules.retrieveAnnotation.output
    output:
        "uORFs/longest_protein_coding_transcripts.gtf"
    conda:
        "envs/uorftools.yaml"
    threads: 1
    shell:
        "mkdir -p uORFs; uORF-Tools/scripts/longest_orf_transcript.py -a {input} -o {output}"

rule normalizedCounts:
    input:
        rules.longestTranscript.output,
        rules.maplink.output
    output:
        "uORFs/sfactors_lprot_FP.csv",
        "uORFs/sfactors_lprot_Total.csv",
        "uORFs/ncounts_lprot_FP.csv",
        "uORFs/ncounts_lprot_Total.csv"
    conda:
        "envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_normalized_counts_longest_protein.R -r -b bam/ -a {input[0]} -s uORFs/sfactors_lprot_FP.csv -n uORFs/ncounts_lprot_FP.csv -t FP;  uORF-Tools/scripts/generate_normalized_counts_longest_protein.R -r -b bam/ -a {input[0]} -s uORFs/sfactors_lprot_Total.csv -n uORFs/ncounts_lprot_Total.csv -t Total")

