import argparse
import subprocess
import os
import sys
import shutil
import warnings
import re

from tqdm import tqdm
from glob import glob

warnings.filterwarnings("ignore")

# ================= Argument Parsing =================
def parse_arguments():
    parser = argparse.ArgumentParser(description="Run the PeGAS pipeline.")
    parser.add_argument("-d", "--data", dest="data", help="Directory containing all the fastq.gz files", required=True)
    parser.add_argument("-o", "--output", dest="output", help="Directory where output files will be saved", required=True)
    parser.add_argument("-c", "--cores", dest="cores", help="The number of cores to use", default=1, type=int)
    # parser.add_argument("--unlock", help="Unlock the working directory", action="store_true")
    parser.add_argument("--overwrite", help="Overwrite the output directory if it exists", action="store_true")
    parser.add_argument("--rerun-pangenome", help="Rerun the pangenome analysis", action="store_true")
    parser.add_argument("--shovill-cpu-cores", dest="shovill_cpu_cores", help="Number of CPU cores to use for Shovill", type=int)
    parser.add_argument("--prokka-cpu-cores", dest="prokka_cpu_cores", help="Number of CPU cores to use for Prokka", type=int)
    parser.add_argument("--roary-cpu-cores", dest="roary_cpu_cores", help="Number of CPU cores to use for Roary", type=int)
    parser.add_argument("--redo-report", help="Debug: redo report", action="store_true")
    return parser.parse_args()

# ================= Utility Functions =================
def list_fastq_files(path):
    """Returns a list of all .fastq.gz files in the specified path relative to base_folder."""
    full_path = os.path.join(path)
    return [f for f in glob(os.path.join(full_path, "*.fastq.gz"))]

def get_core_sample_name(filename):
    """Extracts the core sample name by removing _R1 or _R2 and other suffixes."""
    return os.path.basename(filename).replace("_R1", "").replace("_R2", "").replace(".fastq.gz", "")

def build_fastq_pairs(fastq_files):
    """Pairs R1 and R2 files based on sample names."""
    pairs = {}
    for file in fastq_files:
        sample = get_core_sample_name(file)
        if sample:
            if sample not in pairs:
                pairs[sample] = {}
            if "_R1" in file or "R1" in file:
                pairs[sample]["R1"] = file
            elif "_R2" in file or "R2" in file:
                pairs[sample]["R2"] = file
    return {s: p for s, p in pairs.items() if "R1" in p and "R2" in p}

def copy_files(base_folder, source_files, destination="raw_data"):
    """Copies files from source directory to destination within base_folder without renaming them."""
    dest_path = os.path.join(base_folder, destination)
    os.makedirs(dest_path, exist_ok=True)
    for src in tqdm(source_files, desc="Copying files"):
        dest = os.path.join(dest_path, os.path.basename(src))
        if not os.path.exists(dest) or os.path.getsize(src) != os.path.getsize(dest):
            shutil.copy2(src, dest)
            tqdm.write(f"[pegas] Copied '{src}' to '{dest}'.")

def remove_extra_files(base_folder, destination, valid_files):
    """Removes unwanted files from the destination directory and clears related data for affected samples."""
    dest_path = os.path.join(base_folder, destination)
    valid_basenames = {os.path.basename(f) for f in valid_files}
    valid_samples = list(set([get_core_sample_name(f) for f in valid_files]))

    # Remove extra files in the destination directory
    for file in glob(os.path.join(dest_path, "*")):
        if os.path.basename(file) not in valid_basenames:
            os.remove(file)
            tqdm.write(f"[pegas] Removed '{file}'.")

    # Remove extra sample directories in the results folder
    results_path = os.path.join(base_folder, "results")
    if os.path.exists(results_path):
        for folder in os.listdir(results_path):
            if folder not in valid_samples:
                shutil.rmtree(os.path.join(results_path, folder))
                tqdm.write(f"[pegas] Removed 'results/{folder}' directory.")

    # Remove outdated files in fastqc directory
    fastqc_path = os.path.join(base_folder, "fastqc")
    if os.path.exists(fastqc_path):
        valid_file_names = [os.path.basename(file).replace(".fastq.gz", "") for file in valid_files]
        for file in os.listdir(fastqc_path):
            if os.path.basename(file).replace("_fastqc.html", "").replace("_fastqc.zip", "") not in valid_file_names:
                os.remove(os.path.join(fastqc_path, file))
                tqdm.write(f"[pegas] Removed '{file}'.")

    # Clear obsolete files in pangenome directory
    pangenome_path = os.path.join(base_folder, "pangenome")
    if os.path.exists(pangenome_path):
        for folder in os.listdir(pangenome_path):
            all_files = glob(os.path.join(pangenome_path, folder, "*.gff"))
            for file in all_files:
                if os.path.basename(file).replace(".gff", "") not in valid_samples:
                    tqdm.write(f"[pegas] Removed '{file}' in pangenome.")
                    shutil.rmtree(os.path.join(pangenome_path, folder), ignore_errors=True)
                    
                    # Remove the pangenome flag if needed
                    pangenome_flag_path = os.path.join(base_folder, "flags", ".pangenome")
                    if os.path.exists(pangenome_flag_path):
                        os.remove(pangenome_flag_path)

def main():

    path = os.path.dirname(os.path.realpath(__file__))

    args = parse_arguments()


    data_dir = args.data
    output_dir = args.output
    cores = args.cores
    overwrite = args.overwrite
    rerun_pangenome = args.rerun_pangenome

    # List all FASTQ files in the raw_data_path and raw_data directories
    raw_data_files = list_fastq_files(data_dir)
    # copied_data_files = list_fastq_files(os.path.join(output_dir, "raw_data"))

    # Copy new or modified files from raw_data_path to raw_data
    copy_files(output_dir, raw_data_files, "raw_data")

    # Remove files in raw_data that do not exist in raw_data_path
    remove_extra_files(output_dir, "raw_data", raw_data_files)

    if args.redo_report:
        if os.path.exists(os.path.join(output_dir, "report")):
            shutil.rmtree(os.path.join(output_dir, "report"))
        if os.path.exists(os.path.join(output_dir, "flags", ".report")):
            os.remove(os.path.join(output_dir, "flags", ".report"))

    # Check if the output directory exists
    if os.path.exists(output_dir) and not overwrite:
        tqdm.write("[pegas]Output directory already exists. Use --overwrite to overwrite it or specify a different directory.")
        sys.exit(1)
    else:
        os.makedirs(output_dir, exist_ok=True)

    # Prepare configuration parameters
    config_params = [
        f"raw_data='{data_dir}'",
        f"outdir='{output_dir}'",
        f"install_path='{path}'",
        f"output_dir='{output_dir}'"
    ]

    if rerun_pangenome:
        if os.path.exists(os.path.join(output_dir, "pangenome")):
            shutil.rmtree(os.path.join(output_dir, "pangenome"))
        if os.path.exists(os.path.join(output_dir, "flags", ".pangenome")):
            os.remove(os.path.join(output_dir, "flags", ".pangenome"))

    if args.shovill_cpu_cores:
        config_params.append(f"shovill_cpu_cores={args.shovill_cpu_cores}")
    if args.prokka_cpu_cores:
        config_params.append(f"prokka_cpu_cores={args.prokka_cpu_cores}")
    if args.roary_cpu_cores:
        config_params.append(f"roary_cpu_cores={args.roary_cpu_cores}")

    # Build the Snakemake command
    command = [
        "snakemake",
        "--snakefile", os.path.join(path, "Snakefile"),
        "--directory", output_dir,
        "--cores", str(cores),
        "--rerun-incomplete",
        "--use-conda"
    ]

    unlock_command = [
        "snakemake",
        "--snakefile", os.path.join(path, "Snakefile"),
        "--directory", output_dir,
        "--cores", str(cores),
        "--unlock"
    ]

    # Add the --unlock flag if needed
    # if args.unlock:
    #     command.append("--unlock")

    # Add the --config option and configuration parameters
    command.append("--config")
    command.extend(config_params)

    unlock_command.append("--config")
    unlock_command.extend(config_params)

    tqdm.write("Running PeGAS pipeline with the following command:")
    tqdm.write(" ".join(command))

    # Run the pipeline
    subprocess.run(unlock_command)
    result = subprocess.run(command)
    if result.returncode != 0:
        tqdm.write("Error: Snakemake pipeline failed.")
        sys.exit(result.returncode)
    else:
        tqdm.write("PeGAS pipeline completed successfully.")

if __name__ == "__main__":
    main()