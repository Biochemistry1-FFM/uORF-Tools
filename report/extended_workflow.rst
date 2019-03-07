`uORF-Tools <https://github.com/Biochemistry1-FFM/uORF-Tools>`__ extended workflow performs a analysis of upstream open reading frames (uORF) from ribosome profiling and optionally RNA-seq high troughput sequencing data. The data is processed in the following steps.
Adapter removal is performed with `TrimGalore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`__ and `Cutadapt <http://cutadapt.readthedocs.io>`__, afterwards the reads mapping to rRNA genes are
removed with `SortMeRNA <http://bioinfo.lifl.fr/RNA/sortmerna/>`__. `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__ is used after each of these preprocessing steps. 
The reads are then mapped with `STAR <https://github.com/alexdobin/STAR>`__  
