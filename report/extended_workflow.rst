`uORF-Tools <https://github.com/anibunny12/uORF-Tools>`_ extended workflow performs a analysis of upstream open reading frames (uORF) from ribosome profiling and optionally RNA-seq high troughput sequencing data. The data is processed in the following steps.
Adapter removal is performed with `TrimGalore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`_ and `Cutadapt <http://cutadapt.readthedocs.io>`_, afterwards the reads mapping to rRNA genes are
removed with `SortMeRNA <http://bioinfo.lifl.fr/RNA/sortmerna/>`_. `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ is used after each of these preprocessing steps. 
The reads are then mapped with `STAR <https://github.com/alexdobin/STAR>`_  
Novel ORFs are detected with `RiboTish <https://github.com/zhpn1024/ribotish>`_  
Read counts are normalized with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ .
The normalized reads of detected uORFs together with their main ORF are evaluated regarding 
their ribo_change helping to identify uORFs with a potentially regulatory effect.
