import argparse
import subprocess
import os
import sys
import shutil
import warnings
import json

from tqdm import tqdm
from glob import glob
import multiprocessing as mp

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

try:
    from .gui import launch_gui
    from .common import list_fastq_files, build_fastq_pairs, read_samplesheet, fastq_id, get_core_sample_name, drop_empty_samples
except ImportError:
    from gui import launch_gui
    from common import list_fastq_files, build_fastq_pairs, read_samplesheet, fastq_id, get_core_sample_name, drop_empty_samples

CAP = 16
CONFIG_FILENAME = "pegas_config.json"
TOOL_VERSION_PACKAGES = {
    "fastqc": ["fastqc"],
    "shovill": ["shovill"],
    "amrfinderplus": ["ncbi-amrfinderplus"],
    "abricate": ["abricate"],
    "mlst": ["mlst"],
    "prokka": ["prokka"],
    "quast": ["quast"],
    "checkm2": ["checkm2"],
    "roary": ["roary"],
    "rscript": ["r-base"],
}
CORE_OPTIONS = (
    "shovill_cpu_cores", "prokka_cpu_cores", "amrfinder_cpu_cores",
    "quast_cpu_cores", "checkm2_cpu_cores", "roary_cpu_cores",
)
PROBLEM_STATUSES = {"failed", "skipped", "missing", "no_call"}
OUTPUT_MODES = {
    "new": "new run",
    "resume": "resume: finished steps are reused",
    "overwrite": "start over: previous results will be deleted",
}


def fail(message):
    tqdm.write(f"[pegas] Error: {message}")
    sys.exit(1)


# ---- Argument parsing ----.
def build_parser():
    parser = argparse.ArgumentParser(description="Run the PeGAS pipeline.")
    parser.add_argument("-d", "--data", help="Directory with paired FASTQ files (.fastq.gz, .fq.gz, .fastq, .fq)")
    parser.add_argument("--samplesheet", help="CSV/TSV with columns sample,r1,r2 (alternative to --data)")
    parser.add_argument("-o", "--output", help="Directory for outputs")
    parser.add_argument("-c", "--cores", type=int, help="Total cores to use (default: all)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Delete the previous results in the output directory and start over "
                             "(without it, an existing PeGAS output directory is resumed)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show detected samples and planned resources, then exit")
    parser.add_argument("--setup", action="store_true",
                        help="Create the tool environments, download the databases and check every "
                             "tool starts, then exit (no input needed)")
    parser.add_argument("--shovill-cpu-cores", type=int)
    parser.add_argument("--shovill-ram", type=int, help="RAM (GB) to allocate to Shovill")
    parser.add_argument("--prokka-cpu-cores", type=int)
    parser.add_argument("--amrfinder-cpu-cores", type=int)
    parser.add_argument("--quast-cpu-cores", type=int)
    parser.add_argument("--checkm2-cpu-cores", type=int)
    parser.add_argument("--roary-cpu-cores",  type=int)
    parser.add_argument(
        "--on-tool-failure",
        choices=("stop", "warn"),
        default="warn",
        help="warn (default): mark the failed sample and finish the run for all other samples; "
             "stop: stop the whole run at the first tool failure.",
    )
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
        for option in CORE_OPTIONS[1:]:
            if getattr(args, option) is None:
                setattr(args, option, args.cores)
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
        for option in ("shovill_cpu_cores", "prokka_cpu_cores", "amrfinder_cpu_cores", "quast_cpu_cores"):
            if getattr(args, option) is None:
                setattr(args, option, default_chunk(args.cores))
        for option in ("checkm2_cpu_cores", "roary_cpu_cores"):
            if getattr(args, option) is None:
                setattr(args, option, args.cores)

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

# ---- Inputs and output directory ----.
def resolve_inputs(args):
    """Return (data_dir, output_dir, fastq_files, sample_pairs); exits with a message on unusable input."""
    if args.samplesheet:
        try:
            sample_pairs = read_samplesheet(args.samplesheet)
        except (OSError, ValueError) as exc:
            fail(f"Invalid samplesheet: {exc}")
        data_dir = os.path.dirname(os.path.abspath(os.path.expanduser(args.samplesheet)))
    else:
        if not args.data:
            fail("Provide the input reads with -d/--data or --samplesheet.")
        data_dir = os.path.abspath(os.path.expanduser(args.data))
        if not os.path.isdir(data_dir):
            fail(f"Input directory not found: {data_dir}")
        fastq_files = list_fastq_files(data_dir)
        if not fastq_files:
            fail(f"No FASTQ files (.fastq.gz, .fq.gz, .fastq, .fq) found in: {data_dir}")
        sample_pairs = build_fastq_pairs(fastq_files)
    if not sample_pairs:
        fail("No valid R1/R2 pairs found in the input.")
    if not drop_empty_samples(sample_pairs):
        fail("Every sample has empty read files.")
    fastq_files = sorted(f for reads in sample_pairs.values() for f in reads.values())
    if not args.output:
        fail("Provide the output folder with -o/--output.")
    output_dir = os.path.abspath(os.path.expanduser(args.output))
    return data_dir, output_dir, fastq_files, sample_pairs


def prepare_output_dir(output_dir, overwrite, dry_run=False):
    """Return "new", "resume" or "overwrite" for output_dir and get it ready.

    Only folders that are empty or hold a previous PeGAS run (pegas_config.json)
    are accepted, so a mistyped -o can never delete unrelated data.
    """
    is_pegas_run = os.path.isfile(os.path.join(output_dir, CONFIG_FILENAME))
    if os.path.isdir(output_dir) and os.listdir(output_dir) and not is_pegas_run:
        fail(f"{output_dir} is not empty and is not a PeGAS output directory. Choose a new or empty folder.")
    mode = ("overwrite" if overwrite else "resume") if is_pegas_run else "new"
    if dry_run:
        return mode
    if mode == "overwrite":
        tqdm.write(f"[pegas] --overwrite: deleting previous results in {output_dir}")
        for entry in os.listdir(output_dir):
            if entry == ".snakemake":  # keep Snakemake metadata
                continue
            path = os.path.join(output_dir, entry)
            if os.path.isdir(path) and not os.path.islink(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
    elif mode == "resume":
        tqdm.write(f"[pegas] Resuming the existing run in {output_dir} (use --overwrite to start over).")
    os.makedirs(output_dir, exist_ok=True)
    return mode


def print_plan(args, mode, data_dir, output_dir, sample_pairs):
    say = tqdm.write
    say("[pegas] Dry run: nothing will be executed.")
    say(f"[pegas] Input:  {args.samplesheet or data_dir}")
    say(f"[pegas] Output: {output_dir} ({OUTPUT_MODES[mode]})")
    say(f"[pegas] {len(sample_pairs)} sample(s):")
    for sample, reads in sorted(sample_pairs.items()):
        say(f"[pegas]   {sample}: {os.path.basename(reads['R1'])} + {os.path.basename(reads['R2'])}")
    cores = ", ".join(f"{opt.split('_')[0]} {getattr(args, opt)}" for opt in CORE_OPTIONS)
    say(f"[pegas] Cores: total {args.cores}; {cores}")
    say(f"[pegas] Shovill RAM: {args.shovill_ram or 'auto'} GB; on tool failure: {args.on_tool_failure}")


def print_run_summary(output_dir):
    """Print per-sample problems and QC verdicts from the dataframe outputs."""
    import pandas as pd

    def read_tsv(name):
        try:
            return pd.read_csv(os.path.join(output_dir, "dataframe", name), sep="\t", dtype=str).fillna("")
        except Exception:
            return None

    status, qc = read_tsv("sample_status.tsv"), read_tsv("assembly_qc.tsv")
    if status is None:
        print_report_links(output_dir)
        return
    qc_by_sample = {}
    if qc is not None and "QC_STATUS" in qc.columns:
        qc_by_sample = dict(zip(qc["SAMPLE"], zip(qc["QC_STATUS"], qc["QC_NOTES"])))
    stages = [c for c in status.columns if c not in ("SAMPLE", "NOTES")]
    attention = []
    for _, row in status.iterrows():
        issues = [f"{stage.lower()} {row[stage]}" for stage in stages if row[stage] in PROBLEM_STATUSES]
        qc_status, qc_notes = qc_by_sample.get(row["SAMPLE"], ("", ""))
        if qc_status in ("WARN", "FAIL"):
            issues.append(f"QC {qc_status} ({qc_notes})")
        if issues:
            attention.append(f"{row['SAMPLE']}: {'; '.join(issues)}")
    verdicts = [v for v, _ in qc_by_sample.values()]
    qc_counts = ", ".join(f"{v} {verdicts.count(v)}" for v in ("PASS", "WARN", "FAIL"))
    tqdm.write(f"[pegas] Summary: {len(status)} sample(s); QC {qc_counts}; {len(attention)} need attention.")
    for line in attention[:20]:
        tqdm.write(f"[pegas]   {line}")
    if len(attention) > 20:
        tqdm.write(f"[pegas]   ... and {len(attention) - 20} more (see dataframe/sample_status.tsv)")
    print_report_links(output_dir)


def print_report_links(output_dir):
    """Last lines of a run: file:// links that terminals open on (Ctrl/Cmd+)click."""
    from pathlib import Path
    for report in ("report_r.html", "report.html"):
        path = Path(output_dir, "report", report).resolve()
        if path.is_file():
            tqdm.write(f"[pegas] Open the report: {path.as_uri()}")


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
        "samplesheet": getattr(args, "samplesheet", None),
        "output": output_dir,
        "cores": args.cores,
        "shovill_ram": args.shovill_ram,
        **{option: getattr(args, option, None) for option in CORE_OPTIONS},
        "on_tool_failure": args.on_tool_failure,
        "simple_report": bool(args.simple_report),
        "interactive_report": bool(args.interactive),
    }
    if tool_versions is not None:
        config["tool_versions"] = tool_versions
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    return config_path

def remove_extra_files(base_folder, destination, valid_files, valid_samples):
    """Removes unwanted files from the destination directory and clears related data for affected samples.

    valid_samples are the sample names of this run (the keys of the sample pairs), so a
    samplesheet name like isolate_A is kept even when its files are named A_forward.fq.gz.
    """
    dest_path = os.path.join(base_folder, destination)
    valid_basenames = {os.path.basename(f) for f in valid_files}
    valid_samples = set(valid_samples)

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
        valid_file_names = [fastq_id(file) for file in valid_files]
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

    try:
        from . import pipeline
    except ImportError:
        import pipeline
    if args.setup:
        pipeline.run_setup()
        return

    data_dir, output_dir, raw_data_files, sample_pairs = resolve_inputs(args)
    mode = prepare_output_dir(output_dir, args.overwrite, dry_run=args.dry_run)
    if args.dry_run:
        print_plan(args, mode, data_dir, output_dir, sample_pairs)
        return
    cores = args.cores

    # Remove outputs that do not exist in the input set.
    remove_extra_files(output_dir, "raw_data", raw_data_files, sample_pairs)

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
        f"on_tool_failure='{args.on_tool_failure}'",
        f"envs_dir='{pipeline.ENV_DIR}'",
    ]
    for option in ("shovill_ram",) + CORE_OPTIONS:
        value = getattr(args, option, None)
        if value:
            config_params.append(f"{option}={value}")

    # Build the Snakemake command. mtime-only rerun triggers make a resumed run
    # reuse finished outputs even after a PeGAS upgrade changes rule code.
    command = [
        "snakemake",
        "--snakefile", os.path.join(path, "Snakefile"),
        "--directory", output_dir,
        "--cores", str(cores),
        "--rerun-incomplete",
        "--rerun-triggers", "mtime",
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

    # Tool envs are created (or reused) by the same code as pegas-lite, so both
    # runners and all runs share ~/.pegas/envs.
    pipeline.ensure_conda_available()
    pipeline.ensure_required_envs(args)
    tool_env = pipeline.clean_env()
    subprocess.run(unlock_command, env=tool_env)
    tool_versions = collect_tool_versions([pipeline.ENV_DIR])
    write_run_config(output_dir, data_dir, args, tool_versions=tool_versions)

    # Run the pipeline.
    result = subprocess.run(command, env=tool_env)
    if result.returncode != 0:
        tqdm.write("Error: Snakemake pipeline failed. Per-step logs are in logs/<step>/<sample>.log.")
        sys.exit(result.returncode)
    else:
        tqdm.write("PeGAS pipeline completed successfully.")
        tool_versions = collect_tool_versions([pipeline.ENV_DIR])
        write_run_config(output_dir, data_dir, args, tool_versions=tool_versions)
        print_run_summary(output_dir)

if __name__ == "__main__":
    main()
