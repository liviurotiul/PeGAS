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
    from .main import (
        build_parser,
        finalize_arguments,
        list_fastq_files,
        get_core_sample_name,
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
    from pegas.main import (
        build_parser,
        finalize_arguments,
        list_fastq_files,
        get_core_sample_name,
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
    "quast": "quast_env.yml",
    "checkm2": "checkm2_env.yml",
    "roary": "roary_env.yml",
    "report_r": "report_r_env.yml",
}
ENV_CACHE = {}
CONDA_CHECKED = False
CONDA_PLATFORM = None
# Resolved CheckM2 DIAMOND database path. None = not yet resolved, "" = resolution
# attempted and failed (so we skip CheckM2 without retrying every sample).
CHECKM2_DB = None
CHECKM2_DB_DIR = os.path.join(os.path.expanduser("~"), ".pegas", "checkm2_db")
ROARY_FALLBACK_SPECS = [
    "roary=3.13.0",
    "roary=3.13",
    "roary=3.12",
    "roary",
]


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


def get_conda_platform():
    global CONDA_PLATFORM
    if CONDA_PLATFORM:
        return CONDA_PLATFORM
    conda_subdir = os.environ.get("CONDA_SUBDIR")
    if conda_subdir:
        CONDA_PLATFORM = conda_subdir
        return CONDA_PLATFORM
    try:
        result = subprocess.run(
            ["conda", "info", "--json"],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            info = json.loads(result.stdout)
            platform_name = info.get("platform")
            if platform_name:
                CONDA_PLATFORM = platform_name
                return CONDA_PLATFORM
    except Exception:
        pass
    if sys.platform.startswith("linux"):
        CONDA_PLATFORM = "linux-64"
    elif sys.platform == "darwin":
        CONDA_PLATFORM = "osx-64"
    else:
        CONDA_PLATFORM = sys.platform
    return CONDA_PLATFORM


def get_pegas_version():
    try:
        import importlib.metadata as metadata
        return metadata.version("pegas")
    except Exception:
        pass
    try:
        project_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), "pyproject.toml")
        with open(project_path, "r") as f:
            text = f.read()
        match = re.search(r"^version\s*=\s*['\"]([^'\"]+)['\"]", text, re.MULTILINE)
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


def pin_path_for_yml(yml_path, platform_name=None):
    platform_name = platform_name or get_conda_platform()
    root, _ = os.path.splitext(yml_path)
    candidates = [
        f"{root}.{platform_name}.pin.txt",
        f"{yml_path}.{platform_name}.pin.txt",
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    return None


def explicit_spec_lines_from_text(text):
    lines = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line == "@EXPLICIT":
            continue
        lines.append(line)
    return lines


def read_explicit_spec_lines(path):
    with open(path, "r") as f:
        return explicit_spec_lines_from_text(f.read())


def env_matches_pin(env_path, pin_path):
    if not os.path.isdir(os.path.join(env_path, "conda-meta")):
        return False
    result = subprocess.run(
        ["conda", "list", "-p", env_path, "--explicit", "--md5"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        return False
    return explicit_spec_lines_from_text(result.stdout) == read_explicit_spec_lines(pin_path)


def write_env_meta(meta_path, yml_hash, source_type, source_hash, platform_name, fallback_spec=None):
    meta = {
        "yml_hash": yml_hash,
        "source_type": source_type,
        "source_hash": source_hash,
        "platform": platform_name,
    }
    if fallback_spec is not None:
        meta["fallback_spec"] = fallback_spec
    with open(meta_path, "w") as f:
        json.dump(meta, f)


def env_meta_matches(meta, yml_hash, source_type, source_hash, platform_name):
    if meta.get("yml_hash") != yml_hash:
        return False
    if meta.get("source_type") != source_type:
        return False
    if meta.get("source_hash") != source_hash:
        return False
    if meta.get("platform") != platform_name:
        return False
    return True


def create_roary_env_with_fallbacks(env_path):
    for spec in ROARY_FALLBACK_SPECS:
        if os.path.isdir(env_path):
            shutil.rmtree(env_path, ignore_errors=True)
        tqdm.write(f"[pegas] Retrying roary env with fallback spec '{spec}'.")
        cmd = [
            "conda", "create", "-p", env_path, "-y",
            "-c", "conda-forge",
            "-c", "bioconda",
            "-c", "defaults",
            spec,
        ]
        result = subprocess.run(cmd)
        if result.returncode == 0:
            return spec
    return None


def ensure_env(tool_name):
    ensure_conda_available()
    if tool_name in ENV_CACHE:
        return ENV_CACHE[tool_name]

    env_path, yml_path, meta_path = env_paths(tool_name)
    ensure_dir(os.path.dirname(env_path))

    if not os.path.isfile(yml_path):
        tqdm.write(f"[pegas] Error: missing environment file: {yml_path}")
        sys.exit(1)

    platform_name = get_conda_platform()
    pin_path = pin_path_for_yml(yml_path, platform_name)
    current_hash = hash_file(yml_path)
    pin_hash = hash_file(pin_path) if pin_path else None
    source_type = "pin" if pin_path else "yaml"
    source_hash = pin_hash or current_hash
    if os.path.isdir(env_path) and os.path.isfile(meta_path):
        try:
            with open(meta_path, "r") as f:
                meta = json.load(f)
        except Exception:
            meta = {}
        if env_meta_matches(meta, current_hash, source_type, source_hash, platform_name):
            ENV_CACHE[tool_name] = env_path
            return env_path
        if pin_path and meta.get("yml_hash") == current_hash and env_matches_pin(env_path, pin_path):
            write_env_meta(meta_path, current_hash, "pin", pin_hash, platform_name)
            ENV_CACHE[tool_name] = env_path
            return env_path
        if not pin_path and meta.get("yml_hash") == current_hash and not meta.get("source_type"):
            write_env_meta(meta_path, current_hash, "yaml", current_hash, platform_name)
            ENV_CACHE[tool_name] = env_path
            return env_path

    if os.path.isdir(env_path):
        tqdm.write(f"[pegas] Rebuilding env for {tool_name}.")
        shutil.rmtree(env_path, ignore_errors=True)

    if pin_path:
        tqdm.write(f"[pegas] Creating env for {tool_name} from {platform_name} pin file.")
        cmd = ["conda", "create", "-p", env_path, "--file", pin_path, "-y"]
        result = subprocess.run(cmd)
        if result.returncode == 0:
            write_env_meta(meta_path, current_hash, "pin", pin_hash, platform_name)
            ENV_CACHE[tool_name] = env_path
            return env_path
        tqdm.write(f"[pegas] Warning: pinned env creation failed for {tool_name}; falling back to YAML solve.")
        if os.path.isdir(env_path):
            shutil.rmtree(env_path, ignore_errors=True)

    tqdm.write(f"[pegas] Creating env for {tool_name} from YAML.")
    cmd = ["conda", "env", "create", "-p", env_path, "-f", yml_path, "-y"]
    result = subprocess.run(cmd)
    fallback_spec = None
    if result.returncode != 0:
        if tool_name == "roary":
            fallback_spec = create_roary_env_with_fallbacks(env_path)
        if fallback_spec is None:
            raise subprocess.CalledProcessError(result.returncode, cmd)

    write_env_meta(meta_path, current_hash, "yaml", current_hash, platform_name, fallback_spec=fallback_spec)

    ENV_CACHE[tool_name] = env_path
    return env_path


def _pick_shell():
    """Return a usable shell path for child tools (GNU parallel, prokka, ...)."""
    candidate = os.environ.get("SHELL")
    if candidate and os.path.exists(candidate):
        return candidate
    for fallback in ("/bin/bash", "/usr/bin/bash", "/bin/sh"):
        if os.path.exists(fallback):
            return fallback
    return "/bin/sh"


def tool_env():
    """Environment for tool subprocesses with a guaranteed, exported shell.

    GNU parallel (>=~20250000, e.g. the version prokka shells out to) aborts
    every job with "exec failed: No such file or directory" when it cannot
    resolve a shell from the environment. Under our `conda run -> prokka ->
    perl system() -> /bin/sh -c` chain $SHELL is not exported, so we set both
    SHELL and PARALLEL_SHELL explicitly. See run_prokka().
    """
    env = os.environ.copy()
    shell = _pick_shell()
    env["SHELL"] = shell
    env.setdefault("PARALLEL_SHELL", shell)
    return env


def run_command(cmd, cwd, env_name=None, stdout=None, log_path=None):
    if env_name:
        env_path = ensure_env(env_name)
        cmd = ["conda", "run", "-p", env_path, "--no-capture-output"] + cmd
    if log_path is None:
        result = subprocess.run(cmd, cwd=cwd, stdout=stdout, env=tool_env())
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, cmd)
        return
    # Capture the tool's output to a per-stage log so failures are
    # diagnosable. When the caller redirects stdout to a real output file
    # (e.g. abricate/mlst TSVs), only stderr goes to the log; otherwise the
    # combined stdout+stderr is logged.
    ensure_dir(os.path.dirname(log_path))
    with open(log_path, "w") as log_fh:
        if stdout is None:
            result = subprocess.run(
                cmd, cwd=cwd, stdout=log_fh, stderr=subprocess.STDOUT, env=tool_env()
            )
        else:
            result = subprocess.run(
                cmd, cwd=cwd, stdout=stdout, stderr=log_fh, env=tool_env()
            )
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


def status_dir(work_dir, sample):
    path = os.path.join(work_dir, "results", sample, "status")
    ensure_dir(path)
    return path


def stage_log_path(work_dir, sample, stage):
    # Lives next to the status markers. Ends in .log (not .txt) so it is not
    # picked up by the status-marker parser in build_dataframe.
    return os.path.join(status_dir(work_dir, sample), f"{stage}.log")


def log_tail(log_path, max_lines=25):
    if not log_path or not os.path.isfile(log_path):
        return ""
    try:
        with open(log_path, "r", errors="replace") as fh:
            lines = fh.readlines()
    except OSError:
        return ""
    return "".join(lines[-max_lines:]).rstrip()


def write_stage_status(work_dir, sample, stage, status, message):
    path = os.path.join(status_dir(work_dir, sample), f"{stage}.{status}.txt")
    with open(path, "w") as f:
        f.write(message.rstrip() + "\n")
    return path


def handle_tool_failure(work_dir, sample, stage, exc, failure_policy, log_path=None):
    message = f"{stage} failed for {sample}: {exc}"
    if log_path and os.path.isfile(log_path):
        message += f" (see {os.path.relpath(log_path, work_dir)})"
    write_stage_status(work_dir, sample, stage, "failed", message)
    tail = log_tail(log_path)
    if failure_policy == "warn":
        tqdm.write(f"[pegas] Warning: {message}")
        if tail:
            tqdm.write(f"[pegas] --- {stage} log tail ---\n{tail}")
        return
    tqdm.write(f"[pegas] Error: {message}")
    if tail:
        tqdm.write(f"[pegas] --- {stage} log tail ---\n{tail}")
    raise exc


def valid_assembly(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    try:
        with open(path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    return True
    except OSError:
        return False
    return False


def run_fastqc(fastq_files, work_dir, failure_policy="stop"):
    fastqc_dir = os.path.join(work_dir, "fastqc")
    ensure_dir(fastqc_dir)

    for fastq in tqdm(fastq_files, desc="FastQC"):
        file_id = os.path.basename(fastq).replace(".fastq.gz", "")
        html_path = os.path.join(fastqc_dir, f"{file_id}_fastqc.html")
        zip_path = os.path.join(fastqc_dir, f"{file_id}_fastqc.zip")
        if os.path.exists(html_path) and os.path.exists(zip_path):
            continue
        sample = get_core_sample_name(fastq)
        log_path = stage_log_path(work_dir, sample, "fastqc")
        try:
            run_command(["fastqc", fastq, "-o", fastqc_dir], cwd=work_dir, env_name="fastqc", log_path=log_path)
        except subprocess.CalledProcessError as exc:
            handle_tool_failure(work_dir, sample, "fastqc", exc, failure_policy, log_path=log_path)


def run_shovill(sample_pairs, shovill_cpus, shovill_ram, work_dir, failure_policy="stop"):
    for sample, reads in tqdm(sample_pairs.items(), desc="Shovill"):
        outdir = os.path.join(work_dir, "results", sample, "shovill")
        ensure_dir(outdir)
        assembly = os.path.join(outdir, "contigs.fa")

        if valid_assembly(assembly):
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

        log_path = stage_log_path(work_dir, sample, "shovill")
        try:
            run_command(cmd, cwd=work_dir, env_name="shovill", log_path=log_path)
        except subprocess.CalledProcessError as exc:
            handle_tool_failure(work_dir, sample, "shovill", exc, failure_policy, log_path=log_path)


def run_abricate(sample_names, work_dir, db_name, failure_policy="stop"):
    for sample in tqdm(sample_names, desc=f"Abricate {db_name}"):
        assembly = os.path.join(work_dir, "results", sample, "shovill", "contigs.fa")
        output = os.path.join(work_dir, "results", sample, f"abricate_{db_name}.tsv")
        if os.path.exists(output):
            continue
        if not valid_assembly(assembly):
            write_stage_status(
                work_dir,
                sample,
                f"abricate_{db_name}",
                "skipped",
                f"Abricate {db_name} skipped for {sample}: valid assembly not found.",
            )
            continue
        ensure_dir(os.path.dirname(output))
        log_path = stage_log_path(work_dir, sample, f"abricate_{db_name}")
        try:
            with open(output, "w") as fh:
                run_command(["abricate", "--db", db_name, assembly], cwd=work_dir, env_name="abricate", stdout=fh, log_path=log_path)
        except subprocess.CalledProcessError as exc:
            if os.path.exists(output):
                os.remove(output)
            handle_tool_failure(work_dir, sample, f"abricate_{db_name}", exc, failure_policy, log_path=log_path)


def run_mlst(sample_names, work_dir, failure_policy="stop"):
    for sample in tqdm(sample_names, desc="MLST"):
        assembly = os.path.join(work_dir, "results", sample, "shovill", "contigs.fa")
        output = os.path.join(work_dir, "results", sample, "mlst.tsv")
        if os.path.exists(output):
            continue
        if not valid_assembly(assembly):
            write_stage_status(
                work_dir,
                sample,
                "mlst",
                "skipped",
                f"MLST skipped for {sample}: valid assembly not found.",
            )
            continue
        ensure_dir(os.path.dirname(output))
        log_path = stage_log_path(work_dir, sample, "mlst")
        try:
            with open(output, "w") as fh:
                run_command(["mlst", assembly], cwd=work_dir, env_name="mlst", stdout=fh, log_path=log_path)
        except subprocess.CalledProcessError as exc:
            if os.path.exists(output):
                os.remove(output)
            handle_tool_failure(work_dir, sample, "mlst", exc, failure_policy, log_path=log_path)


def run_prokka(sample_names, prokka_cpus, work_dir, failure_policy="stop"):
    for sample in tqdm(sample_names, desc="Prokka"):
        outdir = os.path.join(work_dir, "results", sample, "prokka")
        ensure_dir(outdir)
        gff = os.path.join(outdir, f"{sample}.gff")
        if os.path.isfile(gff) and os.path.getsize(gff) > 0:
            continue

        assembly = os.path.join(work_dir, "results", sample, "shovill", "contigs.fa")
        if not valid_assembly(assembly):
            write_stage_status(
                work_dir,
                sample,
                "prokka",
                "skipped",
                f"Prokka skipped for {sample}: valid assembly not found.",
            )
            continue
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
        log_path = stage_log_path(work_dir, sample, "prokka")
        try:
            run_command(cmd, cwd=work_dir, env_name="prokka", log_path=log_path)
        except subprocess.CalledProcessError as exc:
            if os.path.isfile(gff) and os.path.getsize(gff) == 0:
                os.remove(gff)
            handle_tool_failure(work_dir, sample, "prokka", exc, failure_policy, log_path=log_path)


def run_quast(sample_names, quast_cpus, work_dir, failure_policy="stop"):
    for sample in tqdm(sample_names, desc="QUAST"):
        outdir = os.path.join(work_dir, "results", sample, "quast")
        report = os.path.join(outdir, "report.tsv")
        if os.path.isfile(report) and os.path.getsize(report) > 0:
            continue

        assembly = os.path.join(work_dir, "results", sample, "shovill", "contigs.fa")
        if not valid_assembly(assembly):
            write_stage_status(
                work_dir,
                sample,
                "quast",
                "skipped",
                f"QUAST skipped for {sample}: valid assembly not found.",
            )
            continue
        ensure_dir(outdir)
        cmd = [
            "quast.py",
            assembly,
            "-o", outdir,
            "--threads", str(quast_cpus),
        ]
        log_path = stage_log_path(work_dir, sample, "quast")
        try:
            run_command(cmd, cwd=work_dir, env_name="quast", log_path=log_path)
        except subprocess.CalledProcessError as exc:
            handle_tool_failure(work_dir, sample, "quast", exc, failure_policy, log_path=log_path)


def _find_checkm2_dmnd(root):
    if not root or not os.path.isdir(root):
        return None
    for current, _dirs, files in os.walk(root):
        for name in files:
            if name.endswith(".dmnd"):
                return os.path.join(current, name)
    return None


def ensure_checkm2_db():
    """Resolve the CheckM2 DIAMOND database, downloading + caching it once.

    Order of preference: $CHECKM2DB if it points at a file, then the cached
    copy under ~/.pegas/checkm2_db, then a one-time `checkm2 database
    --download`. Returns the .dmnd path, or None if it could not be obtained
    (callers then skip CheckM2 rather than failing the run).
    """
    global CHECKM2_DB
    if CHECKM2_DB is not None:
        return CHECKM2_DB or None

    env_db = os.environ.get("CHECKM2DB")
    if env_db and os.path.isfile(env_db):
        CHECKM2_DB = env_db
        return CHECKM2_DB

    cached = _find_checkm2_dmnd(CHECKM2_DB_DIR)
    if cached:
        CHECKM2_DB = cached
        return CHECKM2_DB

    ensure_dir(CHECKM2_DB_DIR)
    env_path = ensure_env("checkm2")
    tqdm.write("[pegas] Downloading CheckM2 database (~3 GB, one-time) ...")
    cmd = [
        "conda", "run", "-p", env_path, "--no-capture-output",
        "checkm2", "database", "--download", "--path", CHECKM2_DB_DIR,
    ]
    result = subprocess.run(cmd, env=tool_env())
    if result.returncode == 0:
        cached = _find_checkm2_dmnd(CHECKM2_DB_DIR)
        if cached:
            CHECKM2_DB = cached
            return CHECKM2_DB
    tqdm.write("[pegas] Warning: CheckM2 database could not be obtained; CheckM2 will be skipped.")
    CHECKM2_DB = ""
    return None


def run_checkm2(sample_names, checkm2_cpus, work_dir, failure_policy="stop"):
    # CheckM2 loads a large DIAMOND database, so all assemblies are scored in a
    # single batch run rather than once per sample.
    targets = []
    for sample in sample_names:
        dest = os.path.join(work_dir, "results", sample, "checkm2", "quality_report.tsv")
        if os.path.isfile(dest) and os.path.getsize(dest) > 0:
            continue
        assembly = os.path.join(work_dir, "results", sample, "shovill", "contigs.fa")
        if not valid_assembly(assembly):
            write_stage_status(
                work_dir,
                sample,
                "checkm2",
                "skipped",
                f"CheckM2 skipped for {sample}: valid assembly not found.",
            )
            continue
        targets.append((sample, assembly))

    if not targets:
        return

    db_path = ensure_checkm2_db()
    if db_path is None:
        for sample, _ in targets:
            write_stage_status(
                work_dir,
                sample,
                "checkm2",
                "skipped",
                f"CheckM2 skipped for {sample}: reference database unavailable.",
            )
        return

    checkm2_dir = os.path.join(work_dir, "checkm2")
    input_dir = os.path.join(checkm2_dir, "input")
    output_dir = os.path.join(checkm2_dir, "output")
    if os.path.isdir(input_dir):
        shutil.rmtree(input_dir, ignore_errors=True)
    ensure_dir(input_dir)

    inputs = []
    for sample, assembly in targets:
        staged = os.path.join(input_dir, f"{sample}.fna")
        try:
            os.symlink(os.path.abspath(assembly), staged)
        except OSError:
            shutil.copy(assembly, staged)
        inputs.append(staged)

    cmd = [
        "checkm2", "predict",
        "--output-directory", output_dir,
        "--threads", str(checkm2_cpus),
        "--database_path", db_path,
        "--force",
        "--input",
    ] + inputs
    log_path = os.path.join(checkm2_dir, "checkm2.log")
    try:
        run_command(cmd, cwd=work_dir, env_name="checkm2", log_path=log_path)
    except subprocess.CalledProcessError as exc:
        tail = log_tail(log_path)
        for sample, _ in targets:
            write_stage_status(
                work_dir, sample, "checkm2", "failed",
                f"checkm2 failed for {sample}: {exc} (see {os.path.relpath(log_path, work_dir)})",
            )
        if failure_policy == "warn":
            tqdm.write(f"[pegas] Warning: checkm2 failed for {len(targets)} sample(s).")
            if tail:
                tqdm.write(f"[pegas] --- checkm2 log tail ---\n{tail}")
            return
        tqdm.write("[pegas] Error: checkm2 failed.")
        if tail:
            tqdm.write(f"[pegas] --- checkm2 log tail ---\n{tail}")
        raise exc

    distribute_checkm2_results(output_dir, targets, work_dir)


def distribute_checkm2_results(output_dir, targets, work_dir):
    report = os.path.join(output_dir, "quality_report.tsv")
    rows_by_sample = {}
    header_line = ""
    if os.path.isfile(report):
        with open(report, "r") as fh:
            lines = fh.read().splitlines()
        if lines:
            header_line = lines[0]
            for line in lines[1:]:
                if not line.strip():
                    continue
                name = line.split("\t", 1)[0]
                rows_by_sample[name] = line

    for sample, _ in targets:
        if sample in rows_by_sample:
            sample_dir = os.path.join(work_dir, "results", sample, "checkm2")
            ensure_dir(sample_dir)
            with open(os.path.join(sample_dir, "quality_report.tsv"), "w") as fh:
                fh.write(header_line + "\n")
                fh.write(rows_by_sample[sample] + "\n")
        else:
            write_stage_status(
                work_dir, sample, "checkm2", "failed",
                f"checkm2 produced no result row for {sample}.",
            )


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
            if not os.path.isfile(src_gff) or os.path.getsize(src_gff) == 0:
                tqdm.write(f"[pegas] Skipping {sample} for Roary: valid Prokka GFF not found.")
                continue
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

        # Remove any stale/partial output before rerunning so Roary does not
        # create a timestamped 'output_*' dir that ends up nested inside it.
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir, ignore_errors=True)

        cmd = [
            "roary",
            "-f", output_dir,
            "-e",
            "-n",
            "-p", str(roary_cpus),
        ] + gff_files
        roary_log = os.path.join(species_dir, "roary.log")
        try:
            run_command(cmd, cwd=work_dir, env_name="roary", log_path=roary_log)
        except subprocess.CalledProcessError:
            tail = log_tail(roary_log)
            if tail:
                tqdm.write(f"[pegas] --- roary log tail ({species}) ---\n{tail}")
            raise

        for item in os.listdir(species_dir):
            if item.startswith("output_") and os.path.isdir(os.path.join(species_dir, item)):
                old_output = os.path.join(species_dir, item)
                subprocess.run(["mv", old_output, output_dir], cwd=work_dir)


def run_report(html_path, data_dir, output_dir, work_dir):
    ensure_dir(os.path.join(work_dir, "report"))
    try:
        from .build_report import build_report
    except ImportError:
        from pegas.build_report import build_report
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
    report_log = os.path.join(work_dir, "report", "report_r.log")
    try:
        run_command(cmd, cwd=work_dir, env_name="report_r", log_path=report_log)
    except subprocess.CalledProcessError:
        tail = log_tail(report_log)
        if tail:
            tqdm.write(f"[pegas] --- report_r log tail ---\n{tail}")
        raise


def write_abricate_database_versions(work_dir):
    ensure_dir(os.path.join(work_dir, "dataframe"))
    output = os.path.join(work_dir, "dataframe", "database_versions.tsv")
    env_path = ensure_env("abricate")
    cmd = ["conda", "run", "-p", env_path, "--no-capture-output", "abricate", "--list"]
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True, check=False)
    if result.returncode == 0 and result.stdout.strip():
        with open(output, "w") as f:
            f.write(result.stdout)
    else:
        with open(output, "w") as f:
            f.write("DATABASE\tSTATUS\tMESSAGE\n")
            message = (result.stderr or result.stdout or "abricate --list failed").strip().replace("\n", " ")
            f.write(f"abricate\tunavailable\t{message}\n")
    return output


def required_envs_for_run(args):
    envs = ["fastqc", "shovill", "abricate", "mlst", "prokka", "quast", "checkm2", "roary"]
    if getattr(args, "simple_report", True):
        envs.append("report_r")
    return envs


def ensure_required_envs(args):
    for env_name in required_envs_for_run(args):
        ensure_env(env_name)


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

    ensure_required_envs(args)

    remove_extra_files(output_dir, "raw_data", fastq_files)
    write_run_config(output_dir, data_dir, args)
    write_sample_manifest(output_dir, fastq_files, sample_pairs)

    work_dir = output_dir
    html_path = os.path.dirname(os.path.realpath(__file__))
    failure_policy = getattr(args, "on_tool_failure", "stop")

    run_fastqc(fastq_files, work_dir, failure_policy=failure_policy)
    run_shovill(sample_pairs, args.shovill_cpu_cores, args.shovill_ram, work_dir, failure_policy=failure_policy)
    sample_names = sorted(sample_pairs.keys())
    run_abricate(sample_names, work_dir, "ncbi", failure_policy=failure_policy)
    run_abricate(sample_names, work_dir, "plasmidfinder", failure_policy=failure_policy)
    run_abricate(sample_names, work_dir, "vfdb", failure_policy=failure_policy)
    run_mlst(sample_names, work_dir, failure_policy=failure_policy)
    run_prokka(sample_names, args.prokka_cpu_cores, work_dir, failure_policy=failure_policy)
    run_quast(sample_names, getattr(args, "quast_cpu_cores", None) or args.cores, work_dir, failure_policy=failure_policy)
    run_checkm2(sample_names, getattr(args, "checkm2_cpu_cores", None) or args.cores, work_dir, failure_policy=failure_policy)

    ensure_dir(os.path.join(work_dir, "dataframe"))
    write_abricate_database_versions(work_dir)
    current_dir = os.getcwd()
    try:
        os.chdir(work_dir)
        build_dataframe()
        run_roary(args.roary_cpu_cores, work_dir)
        tool_versions = collect_tool_versions([ENV_DIR])
        write_run_config(output_dir, data_dir, args, tool_versions=tool_versions)
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
