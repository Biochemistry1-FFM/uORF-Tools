Result table for upstream open reading frames containing eleven columns. 

- uORF_id: uniquely identifies the found uORFs for a transcript by appending an index to the transcript id (e.g ENST00000421495.6.1).

- log2FC_TE_final_uORF: log fold change of overall translational efficiency between conditions for uORF.

- pvalue_final_uORF: significance of differential translation for uORF.

- pvalue.adjust_uORF: significance of differential translation for uORF with multiple testing adjustment.

- transcript_id: transcript identifier of main ORF.

- log2FC_TE_final_CDS: log fold change of overall translational efficiency between conditions for ORF.

- pvalue_final_CDS: significance of differential translation for main ORF.

- pvalue.adjust_CDS: significance of differential translation for main ORF with multiple testing adjustment.

- direction: describes potential regulation of main ORF by uORF by association of their change of translational efficiency.

  - up: translational efficiency of both uORF and main ORF is upregulated.
  
  - down: translational efficency of both uORF and main ORF is downregulated.

  - left: translational efficiency of uORF increases, translational efficiency of main ORF decreases.

  - right: translational efficiency of uORF decreases, translation efficiency of main ORF increases.
  
- ORF_id_gen: genomic coordinates of the potential uORF.

- gene_symbol: gene symbol of main ORF.
