rule riboMerge:
    input:
        expand("ribotish/{sample.condition}-{sample.replicate}-newORFs.tsv", sample=samples.itertuples())
    output:
        "uORFs/Merged_uORF_results.bed",
        "uORFs/Merged_uORF_results.csv"
    conda:
        "../envs/uorftoolspython.yaml"
    threads: 1
    params:
        annotationpath=lambda wildcards: ("NOTSET" if not UORFANNOTATIONPATH else (UORFANNOTATIONPATH))
    shell:
        "if [ -d {params.indexpath} ]; then mkdir -p uORFs; ln -T -s {params.annotationpath} uORFs/Merged_uORF_results.csv; uORF-Tools/scripts/ribo_merge.py --input_csv_filepath {input} --output_bed_filepath uORFs/Merged_uORF_results.bed; else mkdir -p uORFs; uORF-Tools/scripts/ribo_merge.py {input} --min_length 1 --max_length 400 --output_csv_filepath uORFs/Merged_uORF_results.csv --output_bed_filepath uORFs/Merged_uORF_results.bed; fi"

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

rule cdsNormalizedCounts:
    input:
        bam=expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples()),
        annotation=rules.longestTranscript.output,
        sizefactor="uORFs/sfactors_lprot.csv"
    output:
        "uORFs/norm_CDS_reads.csv"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_normalized_counts_CDS.R -b maplink/ -a {input.annotation} -s {input.sizefactor} -t uORF-Tools/samples.tsv -n {output};")

rule uORFNormalizedCounts:
    input:
        bag=expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples()),
        annotation="uORFs/Merged_uORF_results.bed",
        sizefactor="uORFs/sfactors_lprot.csv"
    output:
        "uORFs/norm_uORFs_reads.csv"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_normalized_counts_uORFs.R -b maplink/ -a {input.annotation} -s {input.sizefactor} -t uORF-Tools/samples.tsv -n {output};")

rule cdsxtail:
    input:
        "uORFs/norm_CDS_reads.csv"
    output:
        table=report("uORFs/xtail_cds.csv", caption="../report/xtail_cds_fc.rst", category="CDS"),
        fcplot="uORFs/xtail_cds_fc.pdf",
        rplot="uORFs/xtail_cds_r.pdf"
    conda:
        "../envs/xtail.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/xtail_normalized_counts.R -t uORF-Tools/samples.tsv -r {input} -x {output.table} -f {output.fcplot} -p {output.rplot};")

rule uORFsxtail:
    input:
        "uORFs/norm_uORFs_reads.csv"
    output:
        table=report("uORFs/xtail_uORFs.csv", caption="../report/xtail_uORFs_fc.rst", category="uORFs"),
        tfcplot="uORFs/xtail_uORFs_fc.pdf",
        rplot="uORFs/xtail_uORFs_r.pdf"
    conda:
        "../envs/xtail.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/xtail_normalized_counts.R -t uORF-Tools/samples.tsv -r {input} -x {output.table} -f {output.fcplot} -p {output.rplot};")

rule final_table:
    input:
        xtailuORFs=rules.uORFsxtail.output.table,
	xtailCDS=rules.cdsxtail.output.table,
	annotation="uORFs/Merged_uORF_results.csv"
    output:
        report("uORFs/uORF_regulation.tsv", caption="../report/regulation.rst", category="uORFs")
    conda:
        "../envs/uorftoolspython.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/final_table.py --xtail_cds_file {input.xtailCDS} --xtail_uORF_file {input.xtailuORFs} --uORF_annotation {input.annotation} --output_csv_filepath {output}")
