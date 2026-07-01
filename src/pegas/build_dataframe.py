import json
import os
import re

import pandas as pd
from tqdm import tqdm


RESULT_COLUMNS = [
    "SAMPLE",
    "SPECIES",
    "SUBTYPE",
    "GENE",
    "ACCESSION",
    "PRODUCT",
    "%IDENTITY",
    "%COVERAGE",
    "RESISTANCE",
    "START",
    "END",
    "SEQUENCE",
    "STRAND",
    "PREDICTION_SOURCE",
    "CONTIG_LENGTH",
    "CONTIG_COVERAGE",
    "CONTIG_NUMBER",
]

ABRICATE_COLUMNS = [
    "GENE",
    "ACCESSION",
    "PRODUCT",
    "%IDENTITY",
    "%COVERAGE",
    "RESISTANCE",
    "START",
    "END",
    "SEQUENCE",
    "STRAND",
]

STATUS_COLUMNS = [
    "SAMPLE",
    "FASTQC",
    "SHOVILL",
    "ABRICATE_NCBI",
    "ABRICATE_PLASMIDFINDER",
    "ABRICATE_VFDB",
    "MLST",
    "PROKKA",
    "QUAST",
    "CHECKM2",
    "NOTES",
]

ASSEMBLY_QC_COLUMNS = [
    "SAMPLE",
    "SPECIES",
    "CONTIGS",
    "TOTAL_LENGTH",
    "LARGEST_CONTIG",
    "N50",
    "L50",
    "GC_PERCENT",
    "COMPLETENESS",
    "CONTAMINATION",
]

READ_TOKEN_RE = re.compile(r"R([12])(?=$|[._-])", re.IGNORECASE)


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


def load_manifest_samples():
    manifest_path = "sample_manifest.json"
    if not os.path.isfile(manifest_path):
        return [], {}
    try:
        with open(manifest_path, "r") as f:
            manifest = json.load(f)
    except Exception:
        return [], {}
    sample_pairs = manifest.get("sample_pairs", {})
    if not isinstance(sample_pairs, dict):
        return [], {}
    return sorted(sample_pairs.keys()), sample_pairs


def discover_samples():
    manifest_samples, sample_pairs = load_manifest_samples()
    result_samples = []
    if os.path.isdir("results"):
        result_samples = [
            item for item in os.listdir("results")
            if os.path.isdir(os.path.join("results", item))
        ]
    samples = sorted(set(manifest_samples) | set(result_samples))
    return samples, sample_pairs


def read_stage_markers(sample):
    marker_dir = os.path.join("results", sample, "status")
    markers = {}
    if not os.path.isdir(marker_dir):
        return markers
    for filename in os.listdir(marker_dir):
        if not filename.endswith(".txt"):
            continue
        parts = filename[:-4].split(".")
        if len(parts) != 2:
            continue
        stage, status = parts
        path = os.path.join(marker_dir, filename)
        try:
            with open(path, "r") as f:
                message = f.read().strip()
        except OSError:
            message = ""
        markers[stage] = {"status": status, "message": message}
    return markers


def marker_status(markers, stage):
    entry = markers.get(stage)
    if not entry:
        return None, ""
    return entry.get("status"), entry.get("message", "")


def parse_contigs(sample, markers):
    failed_status, failed_message = marker_status(markers, "shovill")
    contigs_file = os.path.join("results", sample, "shovill", "contigs.fa")
    if failed_status == "failed":
        return 0, 0.0, 0, "failed", failed_message
    if not os.path.isfile(contigs_file):
        return 0, 0.0, 0, "missing", "contigs.fa not found"
    if os.path.getsize(contigs_file) == 0:
        return 0, 0.0, 0, "failed", "contigs.fa is empty"

    try:
        with open(contigs_file, "r") as f:
            contigs_content = f.read()
    except OSError as exc:
        return 0, 0.0, 0, "failed", str(exc)

    contig_entries = contigs_content.strip().split(">")[1:]
    contigs = []  # (length, coverage) pairs kept together so coverage tracks the right contig

    for contig in contig_entries:
        lines = contig.strip().split("\n")
        if not lines or not lines[0]:
            continue
        header = lines[0]
        sequence = "".join(lines[1:])
        len_match = re.search(r"len=(\d+)", header)
        cov_match = re.search(r"cov=([\d.]+)", header)
        length = int(len_match.group(1)) if len_match else len(sequence)
        if length <= 0:
            continue
        coverage = float(cov_match.group(1)) if cov_match else 0.0
        contigs.append((length, coverage))

    if not contigs:
        return 0, 0.0, 0, "failed", "no contigs found in contigs.fa"

    contigs.sort(key=lambda item: item[0], reverse=True)
    half_total = sum(length for length, _ in contigs) / 2
    running = 0
    n50_value, n50_coverage = contigs[-1]
    for length, coverage in contigs:
        running += length
        if running >= half_total:
            n50_value, n50_coverage = length, coverage
            break
    return n50_value, n50_coverage, len(contigs), "ok", ""


def read_mlst(sample, markers):
    failed_status, failed_message = marker_status(markers, "mlst")
    if failed_status in ("failed", "skipped"):
        return "Unknown", "Unknown", failed_status, failed_message

    path = os.path.join("results", sample, "mlst.tsv")
    if not os.path.isfile(path):
        return "Unknown", "Unknown", "missing", "mlst.tsv not found"
    if os.path.getsize(path) == 0:
        return "Unknown", "Unknown", "no_call", "mlst.tsv is empty"
    try:
        mlst = pd.read_csv(path, sep="\t", header=None)
    except Exception as exc:
        return "Unknown", "Unknown", "failed", f"failed to parse mlst.tsv: {exc}"
    if mlst.empty or mlst.shape[1] < 3:
        return "Unknown", "Unknown", "no_call", "MLST returned no species/subtype"
    species = str(mlst.iloc[0, 1]) if pd.notna(mlst.iloc[0, 1]) else "Unknown"
    subtype = str(mlst.iloc[0, 2]) if pd.notna(mlst.iloc[0, 2]) else "Unknown"
    return species, subtype, "ok", ""


def read_quast(sample, markers):
    failed_status, failed_message = marker_status(markers, "quast")
    if failed_status in ("failed", "skipped"):
        return {}, failed_status, failed_message
    path = os.path.join("results", sample, "quast", "report.tsv")
    if not os.path.isfile(path):
        return {}, "missing", "quast report.tsv not found"
    if os.path.getsize(path) == 0:
        return {}, "failed", "quast report.tsv is empty"
    metrics = {}
    try:
        with open(path, "r") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    metrics[parts[0].strip()] = parts[1].strip()
    except OSError as exc:
        return {}, "failed", f"failed to parse quast report.tsv: {exc}"
    out = {
        "CONTIGS": metrics.get("# contigs"),
        "TOTAL_LENGTH": metrics.get("Total length"),
        "LARGEST_CONTIG": metrics.get("Largest contig"),
        "N50": metrics.get("N50"),
        "L50": metrics.get("L50"),
        "GC_PERCENT": metrics.get("GC (%)"),
    }
    return out, "ok", ""


def read_checkm2(sample, markers):
    failed_status, failed_message = marker_status(markers, "checkm2")
    if failed_status in ("failed", "skipped"):
        return {}, failed_status, failed_message
    path = os.path.join("results", sample, "checkm2", "quality_report.tsv")
    if not os.path.isfile(path):
        return {}, "missing", "checkm2 quality_report.tsv not found"
    if os.path.getsize(path) == 0:
        return {}, "failed", "checkm2 quality_report.tsv is empty"
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        return {}, "failed", f"failed to parse checkm2 quality_report.tsv: {exc}"
    if df.empty:
        return {}, "no_call", "checkm2 returned no rows"
    row = df.iloc[0]
    out = {
        "COMPLETENESS": row.get("Completeness"),
        "CONTAMINATION": row.get("Contamination"),
    }
    return out, "ok", ""


def assembly_metrics_from_contigs(sample):
    """Fallback assembly metrics from contigs.fa when QUAST output is absent."""
    path = os.path.join("results", sample, "shovill", "contigs.fa")
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return {}
    lengths = []
    gc = 0
    total = 0
    try:
        with open(path, "r") as f:
            seq_len = 0
            for line in f:
                if line.startswith(">"):
                    if seq_len:
                        lengths.append(seq_len)
                    seq_len = 0
                    continue
                seq = line.strip()
                seq_len += len(seq)
                total += len(seq)
                gc += sum(1 for base in seq.upper() if base in "GC")
            if seq_len:
                lengths.append(seq_len)
    except OSError:
        return {}
    if not lengths:
        return {}
    lengths.sort(reverse=True)
    half = sum(lengths) / 2
    running = 0
    n50 = lengths[-1]
    l50 = len(lengths)
    for idx, length in enumerate(lengths, start=1):
        running += length
        if running >= half:
            n50 = length
            l50 = idx
            break
    gc_pct = (gc / total * 100) if total else None
    return {
        "CONTIGS": len(lengths),
        "TOTAL_LENGTH": total,
        "LARGEST_CONTIG": lengths[0],
        "N50": n50,
        "L50": l50,
        "GC_PERCENT": round(gc_pct, 2) if gc_pct is not None else None,
    }


def read_abricate(sample, db_name, prediction_source, markers):
    stage = f"abricate_{db_name}"
    failed_status, failed_message = marker_status(markers, stage)
    if failed_status in ("failed", "skipped"):
        return pd.DataFrame(columns=ABRICATE_COLUMNS), failed_status, failed_message

    path = os.path.join("results", sample, f"abricate_{db_name}.tsv")
    if not os.path.isfile(path):
        return pd.DataFrame(columns=ABRICATE_COLUMNS), "missing", f"abricate_{db_name}.tsv not found"
    if os.path.getsize(path) == 0:
        return pd.DataFrame(columns=ABRICATE_COLUMNS), "no_hits", ""
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        return pd.DataFrame(columns=ABRICATE_COLUMNS), "failed", f"failed to parse {os.path.basename(path)}: {exc}"

    if df.empty:
        return pd.DataFrame(columns=ABRICATE_COLUMNS), "no_hits", ""
    for column in ABRICATE_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA
    df = df[ABRICATE_COLUMNS].copy()
    df["PREDICTION_SOURCE"] = prediction_source
    return df, "ok", ""


def infer_fastqc_status(sample, sample_pairs, markers):
    failed_status, failed_message = marker_status(markers, "fastqc")
    if failed_status == "failed":
        return "failed", failed_message

    reads = sample_pairs.get(sample, {})
    if not reads:
        return "unknown", "sample not present in manifest"

    missing = []
    for fastq in reads.values():
        file_id = os.path.basename(fastq).replace(".fastq.gz", "")
        html_path = os.path.join("fastqc", f"{file_id}_fastqc.html")
        zip_path = os.path.join("fastqc", f"{file_id}_fastqc.zip")
        if not os.path.exists(html_path) or not os.path.exists(zip_path):
            missing.append(file_id)
    if missing:
        preview = ", ".join(missing[:4])
        more = f" (+{len(missing) - 4} more)" if len(missing) > 4 else ""
        return "missing", f"missing FastQC outputs for {preview}{more}"
    return "ok", ""


def infer_prokka_status(sample, markers):
    failed_status, failed_message = marker_status(markers, "prokka")
    if failed_status in ("failed", "skipped"):
        return failed_status, failed_message
    path = os.path.join("results", sample, "prokka", f"{sample}.gff")
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        return "ok", ""
    if os.path.isfile(path):
        return "failed", "Prokka GFF is empty"
    return "missing", "Prokka GFF not found"


def add_sample_context(df, sample, species, subtype, n50_value, n50_coverage, contig_number):
    if df.empty:
        return df
    df["SAMPLE"] = sample
    df["SPECIES"] = species
    df["SUBTYPE"] = subtype
    df["CONTIG_LENGTH"] = n50_value
    df["CONTIG_COVERAGE"] = n50_coverage
    df["CONTIG_NUMBER"] = contig_number
    return df


def placeholder_row(sample, species, subtype, n50_value, n50_coverage, contig_number):
    row = {column: pd.NA for column in RESULT_COLUMNS}
    row.update({
        "SAMPLE": sample,
        "SPECIES": species,
        "SUBTYPE": subtype,
        "PREDICTION_SOURCE": "None",
        "CONTIG_LENGTH": n50_value,
        "CONTIG_COVERAGE": n50_coverage,
        "CONTIG_NUMBER": contig_number,
    })
    return pd.DataFrame([row], columns=RESULT_COLUMNS)


def finalize_results(df):
    for column in RESULT_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA
    df = df[RESULT_COLUMNS].copy()
    for column in ["%IDENTITY", "%COVERAGE", "START", "END", "CONTIG_LENGTH", "CONTIG_COVERAGE", "CONTIG_NUMBER"]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df["CONTIG_LENGTH"] = df["CONTIG_LENGTH"].fillna(0).astype(int)
    df["CONTIG_NUMBER"] = df["CONTIG_NUMBER"].fillna(0).astype(int)
    df["CONTIG_COVERAGE"] = df["CONTIG_COVERAGE"].fillna(0.0)
    for column in ["SAMPLE", "SPECIES", "SUBTYPE", "GENE", "ACCESSION", "PRODUCT", "RESISTANCE", "SEQUENCE", "STRAND", "PREDICTION_SOURCE"]:
        df[column] = df[column].fillna("").astype(str)
    return df


def build_dataframe():
    os.makedirs("dataframe", exist_ok=True)
    samples, sample_pairs = discover_samples()
    result_frames = []
    status_rows = []
    assembly_qc_rows = []

    for sample in tqdm(samples, desc="Build dataframe"):
        markers = read_stage_markers(sample)
        notes = []

        fastqc_status, fastqc_note = infer_fastqc_status(sample, sample_pairs, markers)
        n50_value, n50_coverage, contig_number, shovill_status, shovill_note = parse_contigs(sample, markers)
        species, subtype, mlst_status, mlst_note = read_mlst(sample, markers)
        prokka_status, prokka_note = infer_prokka_status(sample, markers)

        quast_metrics, quast_status, quast_note = read_quast(sample, markers)
        checkm2_metrics, checkm2_status, checkm2_note = read_checkm2(sample, markers)
        if not any(value not in (None, "") for value in quast_metrics.values()):
            quast_metrics = assembly_metrics_from_contigs(sample)
        assembly_qc_rows.append({
            "SAMPLE": sample,
            "SPECIES": species,
            "CONTIGS": quast_metrics.get("CONTIGS"),
            "TOTAL_LENGTH": quast_metrics.get("TOTAL_LENGTH"),
            "LARGEST_CONTIG": quast_metrics.get("LARGEST_CONTIG"),
            "N50": quast_metrics.get("N50"),
            "L50": quast_metrics.get("L50"),
            "GC_PERCENT": quast_metrics.get("GC_PERCENT"),
            "COMPLETENESS": checkm2_metrics.get("COMPLETENESS"),
            "CONTAMINATION": checkm2_metrics.get("CONTAMINATION"),
        })

        abricate_frames = []
        abricate_statuses = {}
        for db_name, prediction_source, status_column in [
            ("ncbi", "NCBI", "ABRICATE_NCBI"),
            ("plasmidfinder", "PlasmidFinder", "ABRICATE_PLASMIDFINDER"),
            ("vfdb", "VFDB", "ABRICATE_VFDB"),
        ]:
            df_hits, status, note = read_abricate(sample, db_name, prediction_source, markers)
            abricate_statuses[status_column] = status
            if note:
                notes.append(note)
            df_hits = add_sample_context(df_hits, sample, species, subtype, n50_value, n50_coverage, contig_number)
            if not df_hits.empty:
                abricate_frames.append(df_hits)

        for note in [fastqc_note, shovill_note, mlst_note, prokka_note, quast_note, checkm2_note]:
            if note:
                notes.append(note)

        if abricate_frames:
            result_frames.append(pd.concat(abricate_frames, ignore_index=True))
        else:
            result_frames.append(placeholder_row(sample, species, subtype, n50_value, n50_coverage, contig_number))

        status_row = {
            "SAMPLE": sample,
            "FASTQC": fastqc_status,
            "SHOVILL": shovill_status,
            "ABRICATE_NCBI": abricate_statuses["ABRICATE_NCBI"],
            "ABRICATE_PLASMIDFINDER": abricate_statuses["ABRICATE_PLASMIDFINDER"],
            "ABRICATE_VFDB": abricate_statuses["ABRICATE_VFDB"],
            "MLST": mlst_status,
            "PROKKA": prokka_status,
            "QUAST": quast_status,
            "CHECKM2": checkm2_status,
            "NOTES": " | ".join(dict.fromkeys(notes)),
        }
        status_rows.append(status_row)

    if result_frames:
        results = finalize_results(pd.concat(result_frames, ignore_index=True))
    else:
        results = pd.DataFrame(columns=RESULT_COLUMNS)

    status_df = pd.DataFrame(status_rows, columns=STATUS_COLUMNS)

    assembly_qc_df = pd.DataFrame(assembly_qc_rows, columns=ASSEMBLY_QC_COLUMNS)
    for column in ["CONTIGS", "TOTAL_LENGTH", "LARGEST_CONTIG", "N50", "L50", "GC_PERCENT", "COMPLETENESS", "CONTAMINATION"]:
        if column in assembly_qc_df.columns:
            assembly_qc_df[column] = pd.to_numeric(assembly_qc_df[column], errors="coerce")

    results.to_csv("dataframe/results.csv", index=False)
    status_df.to_csv("dataframe/sample_status.tsv", sep="\t", index=False)
    assembly_qc_df.to_csv("dataframe/assembly_qc.tsv", sep="\t", index=False)
