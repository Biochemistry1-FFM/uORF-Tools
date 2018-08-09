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

The retrieval of input files and running the workflow locally and on a server cluster via a queuing system is
working as follows. Create a project directory and change into it:

         mkdir project
         cd project

Retrieve the uORF-Tools from GitHub:

         git clone git@github.com:anibunny12/uORF-Tools.git

The workflow requires a genome sequence (fasta), an annotation file (gtf) and the sequencing results files (fastq).
We recommend retrieving both the genome and the annotation files for mouse and human from [Gencode](https://www.gencodegenes.org/releases/current.html) and for other species from [Ensembl Genomes](http://ensemblgenomes.org/).
Copy the genome and the annotation file into the project folder, decompress then and name them genome.fa and annotation.gtf.

Create a folder fastq and copy your compressed fastq.gz files into the fastq folder.

Please copy the template of the sample sheet and the config file into the project folder.

         cp uORF-Tools/templates/config.yaml .
         cp uORF-Tools/templates/samples.tsv .
       
Customize the config.yaml with the used adapter sequence and optionally with the path to a precomputed
STAR genome index. For correct removal of reads mapping to ribosomal genes please specify the taxonomic group of
the used organism (Eukarya, Bacteria, Archea).
Now edit the sample sheet corresponding to your project, using one line per sequencing result, stating the used
method (RIBO for ribosome profiling, RNA for RNA-seq), the applied condition (e.g. A, B, CTRL, TREAT), the replicate (e.g. 1, 2,..) and the filename. Following is an example:

|method|	condition |replicate|	fastqFile                 |
|------|-----------|---------|--------------------------------|
|RIBO  |	A         |        1|"fastq/FP-treat-1-2.fastq.gz"   |
|RIBO  |	B         |        1|"fastq/FP-ctrl-1-2.fastq.gz"    |
|RNA   |	A         |        1|"fastq/Total-treat-1-2.fastq.gz"|
|RNA   |	B         |        1|"fastq/Total-ctrl-1-2.fastq.gz" |

Now you can start your workflow.

Run Snakemake locally:

         snakemake --use-conda -s Snakefile --configfile config.yaml --directory ${PWD} -j 20 --latency-wait 60 
         

Run Snakemake on the cluster:
Edit cluster.yaml according to your queuing system and cluster hardware. The following example works for Grid Engine:

       snakemake --use-conda -s Snakefile --configfile config.yaml --directory ${PWD} -j 20 --cluster-config uORF-Tools/cluster.yaml --cluster "qsub -N {cluster.jobname} -cwd -q {cluster.qname} -pe {cluster.parallelenvironment} -l {cluster.memory} -o {cluster.logoutputdir} -e {cluster.erroroutputdir} -j {cluster.joinlogs} -M egg@informatik.uni-freiburg.de" --latency-wait 60 

Once the workflow has finished you can request a automatically generated report.html file with the following command:
         
         snakemake --report report.html
