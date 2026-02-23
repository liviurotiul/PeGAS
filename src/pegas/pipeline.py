import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

import pandas as pd
from tqdm import tqdm

try:
    from .build_dataframe import build_dataframe
    from .build_report import build_report
    from .main import (
        build_parser,
        finalize_arguments,
        list_fastq_files,
        build_fastq_pairs,
        write_sample_manifest,
        write_run_config,
        collect_tool_versions,
        remove_extra_files,
    )
except ImportError:
    package_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
    from pegas.build_dataframe import build_dataframe
    from pegas.build_report import build_report
    from pegas.main import (
        build_parser,
        finalize_arguments,
        list_fastq_files,
        build_fastq_pairs,
        write_sample_manifest,
        write_run_config,
        collect_tool_versions,
        remove_extra_files,
    )


def parse_arguments(argv=None, lite_mode=False):
    parser = build_parser()
    args = parser.parse_args(argv)
    return finalize_arguments(args, lite_mode=lite_mode)

ENV_DIR = os.path.join(os.path.expanduser("~"), ".pegas", "envs")
ENV_SPECS = {
    "fastqc": "fastqc_env.yml",
    "shovill": "shovill_env.yml",
    "abricate": "abricate_env.yml",
    "mlst": "mlst_env.yml",
    "prokka": "prokka_env.yml",
    "roary": "roary_env.yml",
    "report_r": "report_r_env.yml",
}
ENV_CACHE = {}
CONDA_CHECKED = False


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def ensure_conda_available():
    global CONDA_CHECKED
    if CONDA_CHECKED:
        return
    if shutil.which("conda") is None:
        tqdm.write("[pegas] Error: conda is not available on PATH. Please install conda and retry.")
        sys.exit(1)
    CONDA_CHECKED = True


def hash_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def get_pegas_version():
    try:
        import importlib.metadata as metadata
        return metadata.version("pegas")
    except Exception:
        pass
    try:
        setup_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "setup.py")
        with open(setup_path, "r") as f:
            text = f.read()
        match = re.search(r"version\\s*=\\s*['\\\"]([^'\\\"]+)['\\\"]", text)
        if match:
            return match.group(1)
    except Exception:
        pass
    return "unknown"


def env_paths(tool_name):
    envs_base = ENV_DIR
    env_path = os.path.join(envs_base, tool_name)
    yml_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "envs", ENV_SPECS[tool_name])
    meta_path = os.path.join(env_path, ".pegas_env.json")
    return env_path, yml_path, meta_path


def ensure_env(tool_name):
    ensure_conda_available()
    if tool_name in ENV_CACHE:
        return ENV_CACHE[tool_name]

    env_path, yml_path, meta_path = env_paths(tool_name)
    ensure_dir(os.path.dirname(env_path))

    if not os.path.isfile(yml_path):
        tqdm.write(f"[pegas] Error: missing environment file: {yml_path}")
        sys.exit(1)

    current_hash = hash_file(yml_path)
    if os.path.isdir(env_path) and os.path.isfile(meta_path):
        try:
            with open(meta_path, "r") as f:
                meta = json.load(f)
        except Exception:
            meta = {}
        if meta.get("yml_hash") == current_hash:
            ENV_CACHE[tool_name] = env_path
            return env_path

    if os.path.isdir(env_path):
        tqdm.write(f"[pegas] Rebuilding env for {tool_name}.")
        shutil.rmtree(env_path, ignore_errors=True)

    tqdm.write(f"[pegas] Creating env for {tool_name}.")
    cmd = ["conda", "env", "create", "-p", env_path, "-f", yml_path, "-y"]
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)

    with open(meta_path, "w") as f:
        json.dump({"yml_hash": current_hash}, f)

    ENV_CACHE[tool_name] = env_path
    return env_path


def run_command(cmd, cwd, env_name=None, stdout=None):
    if env_name:
        env_path = ensure_env(env_name)
        cmd = ["conda", "run", "-p", env_path, "--no-capture-output"] + cmd
    result = subprocess.run(cmd, cwd=cwd, stdout=stdout)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)

def purge_empty_spades(output_dir):
    results_dir = os.path.join(output_dir, "results")
    if not os.path.isdir(results_dir):
        return
    for sample in os.listdir(results_dir):
        sample_dir = os.path.join(results_dir, sample)
        if not os.path.isdir(sample_dir):
            continue
        spades_path = os.path.join(sample_dir, "shovill", "spades.fasta")
        if os.path.isfile(spades_path) and os.path.getsize(spades_path) == 0:
            try:
                shutil.rmtree(sample_dir)
                tqdm.write(f"[pegas] Removed empty Shovill output for {sample} (empty spades.fasta).")
            except Exception as exc:
                tqdm.write(f"[pegas] Warning: failed to remove {sample_dir}: {exc}")


def run_fastqc(fastq_files, work_dir):
    fastqc_dir = os.path.join(work_dir, "fastqc")
    ensure_dir(fastqc_dir)

    for fastq in tqdm(fastq_files, desc="FastQC"):
        file_id = os.path.basename(fastq).replace(".fastq.gz", "")
        html_path = os.path.join(fastqc_dir, f"{file_id}_fastqc.html")
        zip_path = os.path.join(fastqc_dir, f"{file_id}_fastqc.zip")
        if os.path.exists(html_path) and os.path.exists(zip_path):
            continue
        run_command(["fastqc", fastq, "-o", fastqc_dir], cwd=work_dir, env_name="fastqc")


def run_shovill(sample_pairs, shovill_cpus, shovill_ram, work_dir):
    for sample, reads in tqdm(sample_pairs.items(), desc="Shovill"):
        outdir = os.path.join(work_dir, "results", sample, "shovill")
        ensure_dir(outdir)
        assembly = os.path.join(outdir, "contigs.fa")
        gfa = os.path.join(outdir, "contigs.gfa")
        corrections = os.path.join(outdir, "shovill.corrections")
        spades = os.path.join(outdir, "spades.fasta")
        log_path = os.path.join(outdir, "shovill.log")

        if os.path.exists(assembly):
            continue

        cmd = [
            "shovill",
            "--trim",
            "--outdir", outdir,
            "--R1", reads["R1"],
            "--R2", reads["R2"],
            "--force",
            "--cpus", str(shovill_cpus),
        ]
        if shovill_ram:
            cmd.extend(["--ram", str(shovill_ram)])

        try:
            run_command(cmd, cwd=work_dir, env_name="shovill")
        except subprocess.CalledProcessError:
            tqdm.write(f"[pegas] Shovill failed for {sample}. Creating empty outputs.")
            for path in [assembly, gfa, corrections, spades]:
                ensure_dir(os.path.dirname(path))
                open(path, "a").close()
            if not os.path.exists(log_path):
                open(log_path, "a").close()


def run_abricate(sample_names, work_dir, db_name):
    for sample in tqdm(sample_names, desc=f"Abricate {db_name}"):
        assembly = os.path.join(work_dir, "results", sample, "shovill", "contigs.fa")
        output = os.path.join(work_dir, "results", sample, f"abricate_{db_name}.tsv")
        if os.path.exists(output):
            continue
        ensure_dir(os.path.dirname(output))
        with open(output, "w") as fh:
            run_command(["abricate", "--db", db_name, assembly], cwd=work_dir, env_name="abricate", stdout=fh)


def run_mlst(sample_names, work_dir):
    for sample in tqdm(sample_names, desc="MLST"):
        assembly = os.path.join(work_dir, "results", sample, "shovill", "contigs.fa")
        output = os.path.join(work_dir, "results", sample, "mlst.tsv")
        if os.path.exists(output):
            continue
        ensure_dir(os.path.dirname(output))
        with open(output, "w") as fh:
            run_command(["mlst", assembly], cwd=work_dir, env_name="mlst", stdout=fh)


def run_prokka(sample_names, prokka_cpus, work_dir):
    for sample in tqdm(sample_names, desc="Prokka"):
        outdir = os.path.join(work_dir, "results", sample, "prokka")
        ensure_dir(outdir)
        gff = os.path.join(outdir, f"{sample}.gff")
        if os.path.exists(gff):
            continue

        outputs = [
            gff,
            os.path.join(outdir, f"{sample}.err"),
            os.path.join(outdir, f"{sample}.faa"),
            os.path.join(outdir, f"{sample}.ffn"),
            os.path.join(outdir, f"{sample}.fna"),
            os.path.join(outdir, f"{sample}.fsa"),
            os.path.join(outdir, f"{sample}.gbk"),
            os.path.join(outdir, f"{sample}.log"),
            os.path.join(outdir, f"{sample}.sqn"),
            os.path.join(outdir, f"{sample}.tbl"),
            os.path.join(outdir, f"{sample}.tsv"),
            os.path.join(outdir, f"{sample}.txt"),
        ]

        assembly = os.path.join(work_dir, "results", sample, "shovill", "contigs.fa")
        cmd = [
            "prokka",
            "--centre", "X",
            "--compliant",
            assembly,
            "--outdir", outdir,
            "--force",
            "--cpus", str(prokka_cpus),
            "--prefix", sample,
        ]
        try:
            run_command(cmd, cwd=work_dir, env_name="prokka")
        except subprocess.CalledProcessError:
            tqdm.write(f"[pegas] Prokka failed for {sample}. Creating empty outputs.")
            for path in outputs:
                open(path, "a").close()


def link_gff(src_gff, dest_gff):
    if os.path.isdir(dest_gff):
        return False
    if os.path.exists(dest_gff):
        if os.path.realpath(dest_gff) == os.path.realpath(src_gff):
            return True
        if os.path.islink(dest_gff):
            os.unlink(dest_gff)
        else:
            os.remove(dest_gff)
    src_rel = os.path.relpath(src_gff, os.path.dirname(dest_gff))
    os.symlink(src_rel, dest_gff)
    return True


def run_roary(roary_cpus, work_dir):
    df = pd.read_csv(os.path.join(work_dir, "dataframe", "results.csv"), dtype={"SAMPLE": str})
    df_grouped = df.groupby("SPECIES")["SAMPLE"].nunique().reset_index()
    eligible_species = df_grouped.loc[df_grouped["SAMPLE"] > 1, "SPECIES"].tolist()
    if not eligible_species:
        tqdm.write("[pegas] No species eligible for Roary.")
        return

    ensure_dir(os.path.join(work_dir, "pangenome"))

    for species in tqdm(eligible_species, desc="Roary"):
        species_dir = os.path.join(work_dir, "pangenome", species)
        ensure_dir(species_dir)
        output_dir = os.path.join(species_dir, "output")

        samples = df.loc[df["SPECIES"] == species, "SAMPLE"].unique().tolist()
        for sample in samples:
            src_gff = os.path.join(work_dir, "results", sample, "prokka", f"{sample}.gff")
            dest_gff = os.path.join(species_dir, f"{sample}.gff")
            if not link_gff(src_gff, dest_gff):
                tqdm.write(f"[pegas] Skipping {dest_gff}: expected a file, found a directory.")

        gene_presence = os.path.join(output_dir, "gene_presence_absence.csv")
        summary_stats = os.path.join(output_dir, "summary_statistics.txt")
        if os.path.exists(gene_presence) and os.path.exists(summary_stats):
            continue

        gff_files = [os.path.join(species_dir, f) for f in os.listdir(species_dir) if f.endswith(".gff")]
        if len(gff_files) < 2:
            tqdm.write(f"[pegas] Not enough samples for species {species} to run Roary.")
            continue

        cmd = [
            "roary",
            "-f", output_dir,
            "-e",
            "-n",
            "-p", str(roary_cpus),
        ] + gff_files
        run_command(cmd, cwd=work_dir, env_name="roary")

        for item in os.listdir(species_dir):
            if item.startswith("output_") and os.path.isdir(os.path.join(species_dir, item)):
                old_output = os.path.join(species_dir, item)
                subprocess.run(["mv", old_output, output_dir], cwd=work_dir)


def run_report(html_path, data_dir, output_dir, work_dir):
    ensure_dir(os.path.join(work_dir, "report"))
    build_report(html_path, data_dir, output_dir)


def run_r_report(data_dir, output_dir, work_dir, interactive_report=False):
    ensure_dir(os.path.join(work_dir, "report"))
    base_dir = os.path.dirname(os.path.realpath(__file__))
    rmd_template = os.path.join(base_dir, "report_template.Rmd")
    render_script = os.path.join(base_dir, "render_report.R")
    html_output = os.path.join(work_dir, "report", "report_r.html")
    report_html = os.path.join(work_dir, "report", "report.html") if interactive_report else ""
    dataframe_csv = os.path.join(work_dir, "dataframe", "results.csv")
    pegas_version = get_pegas_version()

    cmd = [
        "Rscript",
        render_script,
        "--rmd", rmd_template,
        "--output", html_output,
        "--dataframe_csv", dataframe_csv,
        "--report_html", report_html,
        "--data_dir", data_dir,
        "--output_dir", output_dir,
        "--pegas_version", pegas_version,
        "--pegas_install_dir", base_dir,
    ]
    run_command(cmd, cwd=work_dir, env_name="report_r")


def run_pipeline(args):
    data_dir = os.path.abspath(os.path.expanduser(args.data))
    output_dir = os.path.abspath(os.path.expanduser(args.output))
    ensure_conda_available()

    if getattr(args, "lite_mode", False):
        opts = vars(args)
        ordered_keys = sorted(opts.keys())
        tqdm.write("[pegas] pegas-lite options:")
        for key in ordered_keys:
            tqdm.write(f"[pegas]   {key}={opts[key]}")

    if not os.path.isdir(data_dir):
        tqdm.write(f"[pegas] Input directory not found: {data_dir}")
        sys.exit(1)

    fastq_files = list_fastq_files(data_dir)
    if not fastq_files:
        tqdm.write(f"[pegas] No .fastq.gz files found in: {data_dir}")
        sys.exit(1)

    sample_pairs = build_fastq_pairs(fastq_files)
    if not sample_pairs:
        tqdm.write("[pegas] No valid R1/R2 pairs found in the input folder.")
        sys.exit(1)

    if os.path.exists(output_dir) and not args.overwrite:
        tqdm.write("[pegas]Output directory already exists. Use --overwrite to overwrite it or specify a different directory.")
        sys.exit(1)
    ensure_dir(output_dir)

    if getattr(args, "lite_mode", False):
        purge_empty_spades(output_dir)

    remove_extra_files(output_dir, "raw_data", fastq_files)
    write_run_config(output_dir, data_dir, args)
    write_sample_manifest(output_dir, fastq_files, sample_pairs)

    work_dir = output_dir
    html_path = os.path.dirname(os.path.realpath(__file__))

    run_fastqc(fastq_files, work_dir)
    run_shovill(sample_pairs, args.shovill_cpu_cores, args.shovill_ram, work_dir)
    sample_names = sorted(sample_pairs.keys())
    run_abricate(sample_names, work_dir, "ncbi")
    run_abricate(sample_names, work_dir, "plasmidfinder")
    run_abricate(sample_names, work_dir, "vfdb")
    run_mlst(sample_names, work_dir)
    run_prokka(sample_names, args.prokka_cpu_cores, work_dir)

    ensure_dir(os.path.join(work_dir, "dataframe"))
    current_dir = os.getcwd()
    try:
        os.chdir(work_dir)
        build_dataframe()
        run_roary(args.roary_cpu_cores, work_dir)
        if args.interactive:
            run_report(html_path, data_dir, output_dir, work_dir)
        if args.simple_report:
            run_r_report(data_dir, output_dir, work_dir, interactive_report=args.interactive)
    finally:
        os.chdir(current_dir)

    tool_versions = collect_tool_versions([ENV_DIR])
    write_run_config(output_dir, data_dir, args, tool_versions=tool_versions)


def main(argv=None):
    args = parse_arguments(argv)
    run_pipeline(args)


if __name__ == "__main__":
    main()
