# uORF-Tools 
uORF-Tools are a workflow and a collection of tools for the analysis of 'upstream Open Reading Frames' (short uORFs).
The workflow is based on the workflow management system snakemake and handles installation of all dependencies,
as well as all processings steps. The source code of uORF-Tools is open source and available under the GPL-3 License.
Installation is described below, for usage please refer to the [Userguide](https://uorf-tools.readthedocs.io/en/latest/index.html).

### <u>Installation via bioconda</u>

uORF-Tools requires snakemake (Version=5.4.2), which can be installed with all dependencies via [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Once you have conda installed simply type:

        $ conda create -c conda-forge -c bioconda -n snakemake snakemake==5.4.2

**IMPORTANT!!!** The workflow relies in an older version of conda. To change your conda version type:

        $ conda install -n base conda=4.5.13

        $ conda activate

        $ conda activate snakemake

Create a project directory and change into it:

         $ mkdir project

         $ cd project

Retrieve the uORF-Tools from GitHub:

         $ git clone git@github.com:Biochemistry1-FFM/uORF-Tools.git

Now you can get started. Usage of the workflow is described in the [Userguide](https://uorf-tools.readthedocs.io/en/latest/index.html).
