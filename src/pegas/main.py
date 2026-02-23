import argparse
import subprocess
import os
import sys
import shutil
import warnings
import re
import json

from tqdm import tqdm
from glob import glob
import multiprocessing as mp

warnings.filterwarnings("ignore")

CAP = 16
READ_TOKEN_RE = re.compile(r"R([12])(?=$|[._-])", re.IGNORECASE)
CONFIG_FILENAME = "pegas_config.json"
TOOL_VERSION_PACKAGES = {
    "fastqc": ["fastqc"],
    "shovill": ["shovill"],
    "abricate": ["abricate"],
    "mlst": ["mlst"],
    "prokka": ["prokka"],
    "roary": ["roary"],
    "rscript": ["r-base"],
}

try:
    from .gui import launch_gui
except ImportError:
    from gui import launch_gui

# ---- Argument parsing ----.
def build_parser():
    parser = argparse.ArgumentParser(description="Run the PeGAS pipeline.")
    parser.add_argument("-d", "--data", required=True, help="Directory with fastq.gz files")
    parser.add_argument("-o", "--output", required=True, help="Directory for outputs")
    parser.add_argument("-c", "--cores", type=int, help="Total cores to use (default: all)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite output dir if it exists")
    parser.add_argument("--shovill-cpu-cores", type=int)
    parser.add_argument("--shovill-ram", type=int, help="RAM (GB) to allocate to Shovill")
    parser.add_argument("--prokka-cpu-cores", type=int)
    parser.add_argument("--roary-cpu-cores",  type=int)
    parser.add_argument("--interactive", action="store_true",
                        help="Generate the interactive HTML report (optional)")
    parser.add_argument("--no-r-report", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--html-report", action="store_true",
                        help=argparse.SUPPRESS)
    return parser

def get_total_ram_gb():
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        if pages and page_size:
            return max(1, int(pages * page_size / (1024 ** 3)))
    except Exception:
        return None
    return None

def finalize_arguments(args, lite_mode=False):
    total = mp.cpu_count()
    args.cores = args.cores or total

    # Use one-quarter of total cores, capped, at least one.
    def default_chunk(n):
        return min(CAP, max(1, n // 4))

    if lite_mode:
        ram_gb = args.shovill_ram
        if ram_gb is None:
            ram_gb = get_total_ram_gb()
        if ram_gb is not None and ram_gb <= 0:
            ram_gb = None
        if ram_gb is not None:
            half_ram = max(1, ram_gb // 2)
            tsingle = min(args.cores, CAP, half_ram)
            rsingle = min(ram_gb, max(16, 2 * tsingle))
            args.shovill_cpu_cores = tsingle
            args.shovill_ram = rsingle
        elif args.shovill_cpu_cores is None:
            args.shovill_cpu_cores = min(CAP, args.cores)
        if args.prokka_cpu_cores is None:
            args.prokka_cpu_cores = args.cores
        if args.roary_cpu_cores is None:
            args.roary_cpu_cores = args.cores
    else:
        if args.shovill_cpu_cores is None and args.shovill_ram is None:
            ram_gb = get_total_ram_gb()
            if ram_gb is not None and ram_gb > 0:
                rjob = min(16, max(1, ram_gb // 2))
                j_parallel = max(1, ram_gb // rjob)
                base_threads = max(1, args.cores // j_parallel)
                tjob = max(2, min(8, base_threads))
                if tjob > args.cores:
                    tjob = args.cores
                args.shovill_cpu_cores = tjob
                args.shovill_ram = rjob
        if args.shovill_cpu_cores is None:
            args.shovill_cpu_cores = default_chunk(args.cores)
        if args.prokka_cpu_cores is None:
            args.prokka_cpu_cores = default_chunk(args.cores)
        if args.roary_cpu_cores is None:
            args.roary_cpu_cores = args.cores

    if not hasattr(args, "interactive"):
        args.interactive = False
    args.interactive = bool(args.interactive or getattr(args, "html_report", False))
    args.simple_report = not bool(getattr(args, "no_r_report", False))
    args.lite_mode = lite_mode

    return args

def parse_arguments(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    return finalize_arguments(args, lite_mode=False)

# ---- Utility functions ----.
def list_fastq_files(path):
    """Returns a list of all .fastq.gz files in the specified path."""
    full_path = os.path.abspath(os.path.expanduser(path))
    return [f for f in glob(os.path.join(full_path, "*.fastq.gz"))]

def parse_fastq_read(filename):
    base = os.path.basename(filename)
    if not base.endswith(".fastq.gz"):
        return None
    stem = base[:-9]
    matches = list(READ_TOKEN_RE.finditer(stem))
    if not matches:
        return None
    match = matches[-1]
    sample = stem[:match.start()].rstrip("._-")
    if not sample:
        return None
    return sample, f"R{match.group(1)}"

def get_core_sample_name(filename):
    """Extracts the core sample name by removing _R1 or _R2 and other suffixes."""
    parsed = parse_fastq_read(filename)
    if parsed:
        return parsed[0]
    base = os.path.basename(filename)
    if base.endswith(".fastq.gz"):
        return base[:-9]
    return os.path.splitext(base)[0]

def build_fastq_pairs(fastq_files):
    """Pairs R1 and R2 files based on sample names."""
    pairs = {}
    unmatched = []
    for file in fastq_files:
        parsed = parse_fastq_read(file)
        if not parsed:
            unmatched.append(file)
            continue
        sample, read = parsed
        if sample not in pairs:
            pairs[sample] = {}
        pairs[sample][read] = file

    orphan_samples = sorted([
        sample for sample, reads in pairs.items()
        if "R1" not in reads or "R2" not in reads
    ])
    if unmatched:
        preview = ", ".join(os.path.basename(f) for f in unmatched[:5])
        more = f" (+{len(unmatched) - 5} more)" if len(unmatched) > 5 else ""
        tqdm.write(f"[pegas] Warning: {len(unmatched)} FASTQ files did not match the R1/R2 pattern: {preview}{more}")
    if orphan_samples:
        preview = ", ".join(orphan_samples[:5])
        more = f" (+{len(orphan_samples) - 5} more)" if len(orphan_samples) > 5 else ""
        tqdm.write(f"[pegas] Warning: {len(orphan_samples)} samples missing R1 or R2: {preview}{more}")

    return {s: p for s, p in pairs.items() if "R1" in p and "R2" in p}

def write_sample_manifest(base_folder, fastq_files, sample_pairs):
    """Writes a JSON manifest so Snakemake can avoid re-scanning inputs."""
    manifest_path = os.path.join(base_folder, "sample_manifest.json")
    manifest = {
        "fastq_files": fastq_files,
        "sample_pairs": sample_pairs
    }
    with open(manifest_path, "w") as f:
        json.dump(manifest, f)
    return manifest_path

def _iter_conda_env_paths(prefixes):
    env_paths = []
    for prefix in prefixes:
        if not prefix or not os.path.isdir(prefix):
            continue
        for entry in os.listdir(prefix):
            env_path = os.path.join(prefix, entry)
            if os.path.isdir(env_path) and os.path.isdir(os.path.join(env_path, "conda-meta")):
                env_paths.append(env_path)
    return env_paths


def _conda_list_packages(env_path):
    try:
        result = subprocess.run(
            ["conda", "list", "-p", env_path, "--json"],
            capture_output=True,
            text=True,
            check=False,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    try:
        data = json.loads(result.stdout)
    except Exception:
        return None
    return {pkg.get("name"): pkg.get("version") for pkg in data if isinstance(pkg, dict)}


def collect_tool_versions(conda_prefixes):
    versions = {tool: "unknown" for tool in TOOL_VERSION_PACKAGES}
    if shutil.which("conda") is None:
        return versions
    package_versions = {}
    for env_path in _iter_conda_env_paths(conda_prefixes):
        env_packages = _conda_list_packages(env_path)
        if not env_packages:
            continue
        for name, version in env_packages.items():
            if name and name not in package_versions:
                package_versions[name] = version
    for tool, package_names in TOOL_VERSION_PACKAGES.items():
        for package_name in package_names:
            version = package_versions.get(package_name)
            if version:
                versions[tool] = version
                break
    return versions


def write_run_config(output_dir, data_dir, args, tool_versions=None):
    """Writes the run configuration to a JSON file in the output directory."""
    config_path = os.path.join(output_dir, CONFIG_FILENAME)
    config = {
        "data": data_dir,
        "output": output_dir,
        "cores": args.cores,
        "overwrite": bool(args.overwrite),
        "shovill_cpu_cores": args.shovill_cpu_cores,
        "shovill_ram": args.shovill_ram,
        "prokka_cpu_cores": args.prokka_cpu_cores,
        "roary_cpu_cores": args.roary_cpu_cores,
        "simple_report": bool(args.simple_report),
        "interactive_report": bool(args.interactive),
        "r_report": bool(args.simple_report),
        "html_report": bool(args.interactive),
    }
    if tool_versions is not None:
        config["tool_versions"] = tool_versions
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    return config_path

def remove_extra_files(base_folder, destination, valid_files):
    """Removes unwanted files from the destination directory and clears related data for affected samples."""
    dest_path = os.path.join(base_folder, destination)
    valid_basenames = {os.path.basename(f) for f in valid_files}
    valid_samples = list(set([get_core_sample_name(f) for f in valid_files]))

    # Remove extra files in the destination directory.
    for file in glob(os.path.join(dest_path, "*")):
        if os.path.basename(file) not in valid_basenames:
            os.remove(file)
            tqdm.write(f"[pegas] Removed '{file}'.")

    # Remove extra sample directories in the results folder.
    results_path = os.path.join(base_folder, "results")
    if os.path.exists(results_path):
        for folder in os.listdir(results_path):
            if folder not in valid_samples:
                shutil.rmtree(os.path.join(results_path, folder))
                tqdm.write(f"[pegas] Removed 'results/{folder}' directory.")

    # Remove outdated files in the fastqc directory.
    fastqc_path = os.path.join(base_folder, "fastqc")
    if os.path.exists(fastqc_path):
        valid_file_names = [os.path.basename(file).replace(".fastq.gz", "") for file in valid_files]
        for file in os.listdir(fastqc_path):
            if os.path.basename(file).replace("_fastqc.html", "").replace("_fastqc.zip", "") not in valid_file_names:
                os.remove(os.path.join(fastqc_path, file))
                tqdm.write(f"[pegas] Removed '{file}'.")

    # Clear obsolete files in the pangenome directory.
    pangenome_path = os.path.join(base_folder, "pangenome")
    if os.path.exists(pangenome_path):
        for folder in os.listdir(pangenome_path):
            all_files = glob(os.path.join(pangenome_path, folder, "*.gff"))
            for file in all_files:
                if os.path.basename(file).replace(".gff", "") not in valid_samples:
                    tqdm.write(f"[pegas] Removed '{file}' in pangenome.")
                    shutil.rmtree(os.path.join(pangenome_path, folder), ignore_errors=True)
                    
                    # Remove the pangenome flag if needed.
                    pangenome_flag_path = os.path.join(base_folder, "flags", ".pangenome")
                    if os.path.exists(pangenome_flag_path):
                        os.remove(pangenome_flag_path)

def main():

    path = os.path.dirname(os.path.realpath(__file__))

    if len(sys.argv) == 1:
        gui_args = launch_gui(CONFIG_FILENAME)
        if not gui_args:
            tqdm.write("[pegas] No arguments provided. Exiting.")
            return
        args = parse_arguments(gui_args)
    else:
        args = parse_arguments()


    data_dir = os.path.abspath(os.path.expanduser(args.data))
    output_dir = os.path.abspath(os.path.expanduser(args.output))
    cores = args.cores
    overwrite = args.overwrite

    if not os.path.isdir(data_dir):
        tqdm.write(f"[pegas] Input directory not found: {data_dir}")
        sys.exit(1)

    # List all FASTQ files in the raw_data_path and raw_data directories.
    raw_data_files = list_fastq_files(data_dir)
    if not raw_data_files:
        tqdm.write(f"[pegas] No .fastq.gz files found in: {data_dir}")
        sys.exit(1)

    sample_pairs = build_fastq_pairs(raw_data_files)
    if not sample_pairs:
        tqdm.write("[pegas] No valid R1/R2 pairs found in the input folder.")
        sys.exit(1)

    # Check if the output directory exists.
    if os.path.exists(output_dir) and not overwrite:
        tqdm.write("[pegas]Output directory already exists. Use --overwrite to overwrite it or specify a different directory.")
        sys.exit(1)
    else:
        os.makedirs(output_dir, exist_ok=True)

    # Remove outputs that do not exist in the input set.
    remove_extra_files(output_dir, "raw_data", raw_data_files)

    # Save the run configuration for later reuse.
    write_run_config(output_dir, data_dir, args)

    # Write a manifest of FASTQ files and pairs for Snakemake.
    manifest_path = write_sample_manifest(output_dir, raw_data_files, sample_pairs)
    
    # Prepare configuration parameters.
    config_params = [
        f"raw_data='{data_dir}'",
        f"outdir='{output_dir}'",
        f"install_path='{path}'",
        f"output_dir='{output_dir}'",
        f"sample_manifest='{manifest_path}'",
        f"simple_report={args.simple_report}",
        f"interactive_report={args.interactive}",
        f"r_report={args.simple_report}",
        f"html_report={args.interactive}"
    ]

    if args.shovill_cpu_cores:
        config_params.append(f"shovill_cpu_cores={args.shovill_cpu_cores}")
    if args.shovill_ram:
        config_params.append(f"shovill_ram={args.shovill_ram}")
    if args.prokka_cpu_cores:
        config_params.append(f"prokka_cpu_cores={args.prokka_cpu_cores}")
    if args.roary_cpu_cores:
        config_params.append(f"roary_cpu_cores={args.roary_cpu_cores}")

    # Build the Snakemake command.
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


    # Add the --config option and configuration parameters.
    command.append("--config")
    command.extend(config_params)

    total_ram_gb = get_total_ram_gb()
    if total_ram_gb is not None and total_ram_gb > 0:
        command.extend(["--resources", f"mem_gb={total_ram_gb}"])

    unlock_command.append("--config")
    unlock_command.extend(config_params)

    tqdm.write("Running PeGAS pipeline with the following command:")
    tqdm.write(" ".join(command))

    # Run the pipeline.
    subprocess.run(unlock_command)
    result = subprocess.run(command)
    if result.returncode != 0:
        tqdm.write("Error: Snakemake pipeline failed.")
        sys.exit(result.returncode)
    else:
        tqdm.write("PeGAS pipeline completed successfully.")
        tool_versions = collect_tool_versions([os.path.join(output_dir, ".snakemake", "conda")])
        write_run_config(output_dir, data_dir, args, tool_versions=tool_versions)

if __name__ == "__main__":
    main()
