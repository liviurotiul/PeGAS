import argparse
import subprocess
import os
import sys
from tqdm import tqdm

def main():

    path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser()

    parser.add_argument("--d", "--data", help="The data directory, with all the fastq.gz files", required=True)
    parser.add_argument("--o", "--output", help="The output directory", required=True)
    parser.add_argument("--s", "--samples", help="The path to a text file with the list of samples to be processed, each on a new line", required=False)
    parser.add_argument("--c", "--cores", help="The number of cores to use", required=False, default=1, type=int)
    parser.add_argument("--overwrite", help="Overwrite the output directory if it exists", action="store_true", default=False)
    parser.add_argument("--rerun-pangenome", help="Rerun the pangenome analysis; sometimes pipeline can crash during pangenome analysis; This command will delete all data computed regarding pangenome and start all over", action="store_true", default=False)

    args = parser.parse_args()
    cores = args.c
    overwrite = args.overwrite
    rerun_pangenome = args.rerun_pangenome

    command = f"snakemake --snakefile {path}/Snakefile --cores {cores} --rerun-incomplete --use-conda --config"

    if rerun_pangenome:
        os.system(f"rm -rf {args.o}/pangenome")
        os.system(f"rm -rf {args.o}/flags/.pangenome")

    if args.d:
        command += f" raw_data={args.d}"
    
    if args.o:
        if os.path.exists(args.o) and overwrite == False:
            print("Output directory already exists. Exiting. If you want to overwrite the directory, use the --overwrite argument or remove it manually.")
            sys.exit(1)
        command += f" outdir={args.o}"

    if args.s:
        command += f" samples={args.s}"
    
    command += f" install_path={path}"

    print("Running PeGAS pipeline")

    # Run the pipeline
    subprocess.run(f"{command}", shell=True)

if __name__ == "__main__":
    main()