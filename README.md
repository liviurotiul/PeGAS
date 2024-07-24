

# PeGAS: A Comprehensive Bioinformatic Solution for Pathogenic Bacterial Genomic Analysis

This is PeGAS, a powerful bioinformatic tool designed for the seamless quality control, assembly, and annotation of Illumina paired-end reads specific to pathogenic bacteria. This tool integrates state-of-the-art open-source software to provide a streamlined and efficient workflow, ensuring accurate insights into the genomic makeup of pathogenic microbial strains.

## Key Features:

- **Quality Control:** Utilize industry-standard tools such as FastQC, and Cutadapt to assess and enhance the quality of Illumina paired-end reads.

- **Assembly:** PeGAS uses SPADES for de novo genome assembly, ensuring comprehensive coverage and accurate representation of pathogenic bacterial genomes.

- **Annotation:** We employ abricate for specific gene profiling and prokka for pangenomic analysis

- **Visualization:** PeGAS uses Plotly for interactive data visualisation in the browser
- **Parallel execution:** Using Snakemake as a base, the workflow is mostly parallel allowing for fast execution

## How to Use:


### 1. Prerequisites:

Before using PeGAS, ensure you have the following prerequisites installed:

- **Conda:** If not already installed, follow the instructions [here](https://github.com/conda-forge/miniforge?tab=readme-ov-file) to install conda Miniforge.
- **Mamba:** After the conda installation, open a new terminal (so that the base environement is active and use this command to install Mamba:
	```bash
	(base)user@user:~$ conda install mamba
	```
- **Installation:**

	Clone PeGAS repository to your desired location using Git:

	```bash
	(base)user@user:~/Desired/location$ git clone https://github.com/liviurotiul/PeGAS-snakemake.git
	```
- **Roary:** Install roary:
	```bash
	(base)user@user:sudo apt-get install roary
	```
- **Environment:** Create the conda environemnt you will use for the execution and activate it:
	```bash
	(base)user@user:~/Desired/location$ conda env create -f environment.yaml
	(base)user@user:~/Desired/location$ conda env create -f prokka_env.yaml
	...
	(base)user@user:~/Desired/location$ conda activate PeGAS-snakemake
	```
### 2.  Using PeGAS:
- First set up a folder in a different path to the one where PeGAS was installed
- Copy all your fastq.gz files in the folder with their original names
- Run this command in the terminal:
	```bash
	(base)user@user:~/Desired/location$ snakemake --cores 32 --rerun-incomplete --use-conda --config raw_data=Path/To/Your/Data
	```
- Alternatively you can provide a text file with the name of the samples (one sample per line) by modifying the command as such:
	```bash
	(base)user@user:~/Desired/location$ snakemake --cores 32 --rerun-incomplete --use-conda --config raw_data=Path/To/Your/Data list_of_samples=Path/To/The/File
	```
	**Set the *--cores* parameter according to your machine capabilities.**
### 3.  Visualising the results:
- After the files have been processed and the analysis is completed, you can visualise the results in the path, in the report folder and all the resulting files in the results folder
