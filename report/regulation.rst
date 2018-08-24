Result table for upstream open reading frames containing nine columns. 

- uORF_id: uniqely identifies the found uORFs for a transcript by appending a index to the transcript id (e.g ENST00000421495.6.1).

- log2FC_TE_final_uORF: log fold change of overall translational efficiency between conditions for uORF

- pvalue_final_uORF: significance of differential translation for uORF

- pvalue.adjust_uORF: significance of differential translation for uORF with multiple testing adjustment.

- transcript_id: transcript identifier of main ORF

- log2FC_TE_final_CDS: log fold change of overall translational efficiency between conditions for ORF

- pvalue_final_CDS: signifiance of differential translation for main ORF

- pvalue.adjust_CDS: signigicance of differential translation for main ORF with multiple testing adjustment

- direction: describes potential regulation of main ORF by uORF by association their change of translational efficiency

  - homodirectional: translational efficiency of both ORFs together is up- or downregulated

  - inverse: translational efficiency of uORF increases, translational efficiency of main ORF decreases
  
- ORF_id_gen: genomic coordinates of the potential uORF

- gene_symbol: gene symbol of main ORF
