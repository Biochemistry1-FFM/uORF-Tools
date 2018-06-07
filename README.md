# uORF-Tools
uORF-Tools are a workflow and a collection of tools for the analysis of 'Upstream Open Reading Frames' (short uORFs).
The workflow is based on the workflow management system snakemake and handles installation of all dependencies,
as well as all processings steps. The source code of uORF-Tools is open source and available under the  License.
Installation and basic usage is described below, for more detailed examples please refer to the included [Userguide](Supplement.pdf)

### <u>Installation via bioconda</u>

uORF-Tools can be installed with all dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

         conda create -c bioconda -c conda-forge -n snakemake snakemake 
         
         source activate snakemake

### <u>Basic usage</u>
Using the workflow requires a genome sequence (fasta), an annotation file (gtf) and the sequencing results files (fastq).
We recommend retrieving both the genome and the annotation files from [Ensembl Genomes](http://ensemblgenomes.org/).
The retrieval of input files and running the workflow locally and on a server cluster via a queuing system is
demonstrated using an example from [ENA](https://www.ebi.ac.uk/ena)

