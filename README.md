
  
  

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

  

<img  src="SunburstDemo1.gif"  width="100%">

  

## How to Use

  
  

### 1. Prerequisites

  

- **CONDA**

  

**IMPORTANT**: PeGAS uses **snakemake** as a base which needs **mamba**, a reimplementation of conda in C++. As far as we know, mamba can **only** be installed alongside miniconda or miniforge. The installation alongside anaconda is not yet possible.

  

To find out if you have miniconda or anaconda you can open a terminal and type:

  

```
(base)user@user:~$ conda info
```

  

You should be able to see either miniconda or anaconda.

  

If you **don't have conda yet installed** we recommend conda [[Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file)](https://conda-forge.org/download/)

  

- **Mamba:** After the conda installation, open a new terminal (so that the base environement is active and use this command to install Mamba:

```bash

(base)user@user:~$ conda install conda-forge::mamba

```

  

### 2. Installing PeGAS

  

- PeGAS is easily installable with conda; just open a terminal and run this command:

  
  
```
(base)user@user:~$ conda install bioconda::pegas
```
or
```
(base)user@user:~$ conda create -n pegas -c bioconda -c conda-forge
```

### 3. Using PeGAS

  

- First set up a folder in a different path to the one where PeGAS was installed

- Copy all your fastq.gz files in the folder with their original names

  

- Pegas should then work almost out of the box just by activating the **pegas** environment and running:

  

```bash

usage: pegas [-h] -d DATA -o OUTPUT [-c CORES] [--overwrite] [--shovill-cpu-cores SHOVILL_CPU_CORES] [--prokka-cpu-cores PROKKA_CPU_CORES] [--roary-cpu-cores ROARY_CPU_CORES] [--gc GC]

Run the PeGAS pipeline.

options:
  -h, --help            show this help message and exit
  -d DATA, --data DATA  Directory containing all the fastq.gz files
  -o OUTPUT, --output OUTPUT
                        Directory where output files will be saved
  -c CORES, --cores CORES
                        The number of cores to use
  --overwrite           Overwrite the output directory if it exists
  --shovill-cpu-cores SHOVILL_CPU_CORES
                        Number of CPU cores to use for Shovill
  --prokka-cpu-cores PROKKA_CPU_CORES
                        Number of CPU cores to use for Prokka
  --roary-cpu-cores ROARY_CPU_CORES
                        Number of CPU cores to use for Roary
  --gc GC               Provide a custom json file path for GC content limits for each species

```
A model for the json file can be found [here](https://github.com/liviurotiul/PeGAS/blob/main/src/pegas/gc_content.json)
### 3. Visualising the results

- After the files have been processed and the analysis is completed, you can visualise the results in the path, in the report folder and all the resulting files in the results folder

### PeGAS works out of the box on linux systems. For Windows and MAC users we recommend using virtualization software like Docker or WSL! 
