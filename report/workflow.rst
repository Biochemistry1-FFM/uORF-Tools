`uORF-Tools <https://github.com/Biochemistry1-FFM/uORF-Tools>`__ performs an analysis of upstream open reading frames (uORFs).
Novel ORFs are detected with `RiboTISH <https://github.com/zhpn1024/ribotish>`__ . 
Read counts are normalized with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`__ .
The normalized reads of detected uORFs together with their main ORFs are evaluated regarding 
their ribo_change helping to identify uORFs with a potentially regulatory effect. 
