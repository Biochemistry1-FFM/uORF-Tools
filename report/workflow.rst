`uORF-Tools <https://github.com/Biochemistry1-FFM/uORF-Tools>`_ performs a analysis of upstream open reading frames (uORF).
Novel ORFs are detected with `RiboTish <https://github.com/zhpn1024/ribotish>`_ . 
Read counts are normalized with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ .
The normalized reads of detected uORFs together with their main ORF are evaluated regarding 
their ribo_change helping to identify uORFs with a potentially regulatory effect. 
