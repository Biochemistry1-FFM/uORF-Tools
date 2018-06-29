`uORF-Tools <https://github.com/anibunny12/uORF-Tools>` performs a analysis of upstream open reading frames (uORF) from ribosome profiling and
RNA-seq high troughput sequencing data. The data is processed in the following steps.
Adapter removal is performed with `TrimGalore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>` and `Cutadapt <http://cutadapt.readthedocs.io>` , afterwards the reads mapping to rRNA genes are
removed with `SortMeRNA <http://bioinfo.lifl.fr/RNA/sortmerna/>` . `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>` is used after each of these preprocessing steps. 
The reads are then mapped with `STAR <https://github.com/alexdobin/STAR>`. New ORFs not contained in the
provided annotation are detected with `RiboTaper <https://ohlerlab.mdc-berlin.de/software/RiboTaper_126/>`. The mapped reads are visualized as
wig files. Diffential expression of ORFs is evaluated with `Xtail <https://github.com/xryanglab/xtail>`. The regulatory
impact of uORFs on their associated ORFs is detected via novel scripts included in uORF-Tools.
