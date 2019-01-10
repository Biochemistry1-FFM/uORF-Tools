rule riboMerge:
    input:
        expand("ribotish/{sample.condition}-{sample.replicate}-newORFs.tsv_all.txt", sample=samples.itertuples())
    output:
        bed="uORFs/merged_uORFs.bed",
        csv="uORFs/merged_uORFs.csv"
    conda:
        "../envs/uorftoolspython.yaml"
    threads: 1
    params:
        annotationpath=lambda wildcards: ("NOTSET" if not UORFANNOTATIONPATH else (UORFANNOTATIONPATH))
    shell:
        "if [ -e {params.annotationpath} ]; then export APATH=`readlink -f {params.annotationpath}`; mkdir -p uORFs; ln -T -s $APATH {output.csv}; uORF-Tools/scripts/ribo_convert.py --input_csv_filepath $APATH --output_bed_filepath {output.bed}; else mkdir -p uORFs; uORF-Tools/scripts/ribo_merge.py {input} --min_length 1 --max_length 400 --output_csv_filepath {output.csv} --output_bed_filepath {output.bed}; fi"

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
        longestTranscript=rules.longestTranscript.output,
        bams=expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples())
    output:
        "uORFs/sfactors_lprot.csv"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_size_factors.R -t uORF-Tools/samples.tsv -b maplink/ -a {input.longestTranscript} -s uORFs/sfactors_lprot.csv;")

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
        annotation="uORFs/merged_uORFs.bed",
        sizefactor="uORFs/sfactors_lprot.csv"
    output:
        "uORFs/norm_uORFs_reads.csv"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_normalized_counts_uORFs.R -b maplink/ -a {input.annotation} -s {input.sizefactor} -t uORF-Tools/samples.tsv -n {output};")

rule cdsRiboCounts:
    input:
        bam=expand("maplink/RIBO/{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples()),
        annotation=rules.longestTranscript.output,
        sizefactor="uORFs/sfactors_lprot.csv"
    output:
        norm="uORFs/ribo_norm_CDS_reads.csv",
        raw="uORFs/ribo_raw_CDS_reads.csv"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_ribo_counts_CDS.R -b maplink/RIBO/ -a {input.annotation} -s {input.sizefactor} -t uORF-Tools/samples.tsv -n {output.norm} -r {output.raw};")

rule uORFRiboCounts:
    input:
        bam=expand("maplink/RIBO/{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples()),
        annotation="uORFs/merged_uORFs.bed",
        sizefactor="uORFs/sfactors_lprot.csv"
    output:
        norm="uORFs/ribo_norm_uORFs_reads.csv",
        raw="uORFs/ribo_raw_uORFs_reads.csv"
    conda:
        "../envs/uorftools.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/generate_ribo_counts_uORFs.R -b maplink/RIBO/ -a {input.annotation} -s {input.sizefactor} -t uORF-Tools/samples.tsv -n {output.norm} -r {output.raw};")

rule riboChanges:
    input:
        uorf="uORFs/ribo_norm_uORFs_reads.csv",
        orf="uORFs/ribo_norm_CDS_reads.csv"
    output:
        frac="uORFs/ribo_change.csv"
    conda:
        "../envs/uorftoolspython.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/ribo_changes.py --uORF_reads {input.uorf} --ORF_read {input.orf} --changes_output {output.frac};")


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
        fcplot="uORFs/xtail_uORFs_fc.pdf",
        rplot="uORFs/xtail_uORFs_r.pdf"
    conda:
        "../envs/xtail.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/xtail_normalized_counts.R -t uORF-Tools/samples.tsv -r {input} -x {output.table} -f {output.fcplot} -p {output.rplot};")

rule final_table:
    input:
        annotation="uORFs/merged_uORFs.csv",
        uORFreads="uORFs/ribo_norm_uORFs_reads.csv",
        cdsreads="uORFs/ribo_norm_CDS_reads.csv"
    output:
        report("uORFs/uORFs_regulation.tsv", caption="../report/regulation.rst", category="uORFs")
    conda:
        "../envs/uorftoolspython.yaml"
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/final_table.py --uORF_reads {input.uORFreads} --ORF_reads {input.cdsreads} --uORF_annotation {input.annotation} --output_csv_filepath {output}")

rule processing_summary:
    input:
        "uORFs/merged_uORFs.csv"
    output:
        report("uORFs/processing_summary.tsv", caption="../report/summary.rst", category="uORFs")
    conda:
        "../envs/uorftoolspython.yaml"
    params:
       cwd=os.getcwd()
    threads: 1
    shell: ("mkdir -p uORFs; uORF-Tools/scripts/summary_results.py --sample_tsv uORF-Tools/samples.tsv --project_folder {params.cwd} --output_tsv_filepath {output}")
