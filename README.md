# uORF-Tools
uORF-Tools are a workflow and a collection of tools for the analysis of 'Upstream Open Reading Frames' (short uORFs).
The workflow is based on the workflow management system snakemake and handles installation of all dependencies,
as well as all processings steps. The source code of uORF-Tools is open source and available under the  License.
Installation and basic usage is described below, for more detailed examples please refer to the included [Userguide](Supplement.pdf)

### <u>Installation via bioconda</u>

uORF-Tools can be installed with all dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

         conda create -c bioconda -c conda-forge -n uORF-Tools snakemake 
         
         source activate uORF-Tools

### <u>Basic usage</u>
Using the workflow requires a genome sequence (fasta), an annotation file (gtf) and the sequencing results files (fastq).
We recommend retrieving both the genome and the annotation files from [Ensembl Genomes](http://ensemblgenomes.org/).
The retrieval of input files and running the workflow locally and on a server cluster via a queuing system is
demonstrated using an example from [ENA](https://www.ebi.ac.uk/ena)

Retrieve uORF-Tools:

         git clone git@github.com:anibunny12/uORF-Tools.git
         
Copy the genome and the annotation file into the uORF-Tools folder and name them genome.fa and annotation.gtf.

Create a folder fastq/ and copy your fastq files into the folder.

Edit sample.tsv and enter the names of your fastq-files. Trim Galore will try to auto-detect the used adapter-sequences,
if known add your adapter sequence to config.yaml.

Run Snakemake locally:

         snakemake --use-conda -s Snakefile --configfile config.yaml --directory ${PWD} -j 20 --latency-wait 60 
         

Run Snakemake on the cluster:
Edit cluster.yaml according to your queuing system and cluster hardware. The following example works for Grid Engine:

       snakemake --use-conda -s Snakefile --configfile config.yaml --directory ${PWD} -j 20 --cluster-config uORF-Tools/cluster.yaml --cluster "qsub -N {cluster.jobname} -cwd -q {cluster.qname} -pe {cluster.parallelenvironment} -l {cluster.memory} -o {cluster.logoutputdir} -e {cluster.erroroutputdir} -j {cluster.joinlogs} -M egg@informatik.uni-freiburg.de" --latency-wait 60 



         

