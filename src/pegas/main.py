import argparse
import subprocess
import os
import sys
from tqdm import tqdm

def main():

    path = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser()

    parser.add_argument("--d", "--data", help="The data directory", required=True)
    parser.add_argument("--o", "--output", help="The output directory", required=True)
    parser.add_argument("--s", "--samples", help="The path to a text file with the list of samples to be processed", required=False)

    args = parser.parse_args()

    command = f"snakemake --snakefile {path}/Snakefile --cores 32 --rerun-incomplete --use-conda --config"

    if args.d:
        command += f" raw_data={args.d}"
    
    if args.o:
        command += f" outdir={args.o}"
    
    if args.s:
        command += f" samples={args.s}"
    
    command += f" install_path={path}"

    print("Running PeGAS pipeline")

    # Run the pipeline
    subprocess.run(f"{command}", shell=True)

if __name__ == "__main__":
    main()