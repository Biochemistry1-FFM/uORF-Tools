# uORF-Tools [![ReadTheDocs](https://readthedocs.org/projects/uorf-tools/badge/?version=latest)](https://uorf-tools.readthedocs.io/en/latest/index.html) [![GitHub](https://img.shields.io/github/tag/Biochemistry1-FFM/uORF-Tools.svg)](https://github.com/Biochemistry1-FFM/uORF-Tools)
uORF-Tools are a workflow and a collection of tools for the analysis of 'Upstream Open Reading Frames' (short uORFs).
The workflow is based on the workflow management system snakemake and handles installation of all dependencies,
as well as all processings steps. The source code of uORF-Tools is open source and available under the GPL-3 License.
Installation is described below, for usage examples and detailed a detailed manual please refer to the [Userguide](https://uorf-tools.readthedocs.io/en/latest/index.html).

### <u>Installation via bioconda</u>

uORF-Tools require snakemake (Version >=5.3.1), which can be installed with all dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

         conda create -c conda-forge -c bioconda -n snakemake snakemake
         
         source activate snakemake

Create a project directory and change into it:

         mkdir project
         cd project

Retrieve the uORF-Tools from GitHub:

         git clone git@github.com:Biochemistry1-FFM/uORF-Tools.git

Using the workflow is described in the [Userguide](https://uorf-tools.readthedocs.io/en/latest/index.html).
