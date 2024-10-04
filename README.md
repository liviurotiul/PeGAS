

# PeGAS: A Comprehensive Bioinformatic Solution for Pathogenic Bacterial Genomic Analysis

This is PeGAS, a powerful bioinformatic tool designed for the seamless quality control, assembly, and annotation of Illumina paired-end reads specific to pathogenic bacteria. This tool integrates state-of-the-art open-source software to provide a streamlined and efficient workflow, ensuring accurate insights into the genomic makeup of pathogenic microbial strains.

<br/>

## Key Features


- **Quality Control:** Utilize industry-standard tools such as FastQC, and Cutadapt to assess and enhance the quality of Illumina paired-end reads.

- **Assembly:** PeGAS uses SPADES for de novo genome assembly, ensuring comprehensive coverage and accurate representation of pathogenic bacterial genomes.

- **Annotation:** We employ abricate for specific gene profiling and prokka for pangenomic analysis

- **Visualization:** PeGAS uses Plotly for interactive data visualisation in the browser

- **Parallel execution:** Using Snakemake as a base, the workflow is mostly parallel allowing for fast execution

![Alt text](Features.png)

<img src="SunburstDemo1.gif" width="100%">

## How to Use


### 1. Prerequisites

- **CONDA**

	**IMPORTANT**: PeGAS uses **snakemake** as a base which needs **mamba**, a reimplementation of conda in C++. As far as we know, mamba can **only** be installed alongside miniconda or miniforge. The installation alongside anaconda is not yet possible.

	To find out if you have miniconda or anaconda you can open a terminal and type:

		conda info

	You should be able to see either miniconda or anaconda a lot.

	If you **don't have conda yet installed** we recommend conda [Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file)

- **Mamba:** After the conda installation, open a new terminal (so that the base environement is active and use this command to install Mamba:
	```bash
	(base)user@user:~$ conda install mamba
	```

### 2. Installing PeGAS



### 3.  Using PeGAS

- First set up a folder in a different path to the one where PeGAS was installed
- Copy all your fastq.gz files in the folder with their original names

- Pegas should then work almost out of the box just by activating the **pegas** environment and running:

	```bash
	options:
	-h, --help            show this help message and exit
	--d D, --data D       The data directory, with all the fastq.gz files
	--o O, --output O     The output directory
	--s S, --samples S    The path to a text file with the list of samples to be processed, each on a new line
	--c C, --cores C      The number of cores to use
	--overwrite           Overwrite the output directory if it exists
	--rerun-pangenome     Rerun the pangenome analysis; sometimes pipeline can crash during pangenome analysis; This command will delete all data computed regarding pangenome and start all over
	--shovill-cpu-cores SHOVILL_CPU_CORES
							Number of CPU cores to use for shovill
	--prokka-cpu-cores PROKKA_CPU_CORES
							Number of CPU cores to use for prokka
	--roary-cpu-cores ROARY_CPU_CORES
							Number of CPU cores to use for roary
	```
### 3.  Visualising the results
- After the files have been processed and the analysis is completed, you can visualise the results in the path, in the report folder and all the resulting files in the results folder
