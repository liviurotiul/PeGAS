"""Helpers shared by the CLI, the lite runner, the Snakefile and build_dataframe."""
import csv
import os
import re
import sys
from glob import glob

try:
    from .utils import species_dict
except ImportError:
    from utils import species_dict

# Order matters: the longest suffix must be tried first. FastQC strips the same
# suffixes when naming its reports, so fastq_id() matches FastQC output names.
FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
# Mate token: "R1"/"R2" before a separator or the end ("S1_R1_001", "sampleR2"),
# or a bare "1"/"2" right after a separator ("SRR123_1"). The last match wins.
READ_TOKEN_RE = re.compile(r"(?:R|(?<=[._-]))([12])(?=$|[._-])", re.IGNORECASE)
SAMPLE_NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
# MLST species values that must never be grouped together for the pangenome.
UNASSIGNED_SPECIES = {"", "-", "nan", "Unknown", "Unassigned"}

# AMRFinderPlus --organism values (enables point-mutation screening), keyed by the
# species names in utils.species_dict. Neisseria is left out on purpose: the MLST
# "neisseria" scheme covers several species with different resistance mutations.
AMRFINDER_ORGANISMS = {name: name.replace(" ", "_") for name in [
    "Acinetobacter baumannii", "Burkholderia cepacia", "Burkholderia pseudomallei",
    "Citrobacter freundii", "Clostridioides difficile", "Enterobacter cloacae",
    "Enterococcus faecalis", "Enterococcus faecium", "Klebsiella oxytoca",
    "Klebsiella pneumoniae", "Pseudomonas aeruginosa", "Staphylococcus aureus",
    "Staphylococcus pseudintermedius", "Streptococcus agalactiae",
    "Streptococcus pneumoniae", "Streptococcus pyogenes", "Vibrio cholerae",
    "Vibrio parahaemolyticus", "Vibrio vulnificus",
]}
AMRFINDER_ORGANISMS.update({
    "Escherichia coli": "Escherichia",
    "Salmonella enterica": "Salmonella",
    "Campylobacter jejuni": "Campylobacter",
})


def warn(message):
    print(f"[pegas] Warning: {message}", file=sys.stderr)


def fastq_id(path):
    """File name without its FASTQ extension, or None if it is not a FASTQ file."""
    base = os.path.basename(path)
    for ext in FASTQ_EXTENSIONS:
        if base.endswith(ext):
            return base[: -len(ext)]
    return None


def list_fastq_files(path):
    full_path = os.path.abspath(os.path.expanduser(path))
    return sorted(f for f in glob(os.path.join(full_path, "*")) if os.path.isfile(f) and fastq_id(f))


def parse_fastq_read(filename):
    """Return (sample, "R1"|"R2") for a FASTQ file name, or None."""
    stem = fastq_id(filename)
    if not stem:
        return None
    matches = list(READ_TOKEN_RE.finditer(stem))
    if not matches:
        return None
    match = matches[-1]
    sample = stem[:match.start()].rstrip("._-")
    return (sample, f"R{match.group(1)}") if sample else None


def get_core_sample_name(filename):
    parsed = parse_fastq_read(filename)
    if parsed:
        return parsed[0]
    return fastq_id(filename) or os.path.splitext(os.path.basename(filename))[0]


def _preview(items, limit=5):
    more = f" (+{len(items) - limit} more)" if len(items) > limit else ""
    return ", ".join(items[:limit]) + more


def build_fastq_pairs(fastq_files):
    """Pair R1/R2 files by sample name; unpaired or unrecognised files are reported and dropped."""
    pairs = {}
    unmatched = []
    duplicates = []
    for file in fastq_files:
        parsed = parse_fastq_read(file)
        if not parsed:
            unmatched.append(os.path.basename(file))
            continue
        sample, read = parsed
        if read in pairs.setdefault(sample, {}):
            duplicates.append(os.path.basename(file))
            continue
        pairs[sample][read] = file
    orphans = sorted(s for s, reads in pairs.items() if len(reads) != 2)
    if unmatched:
        warn(f"{len(unmatched)} FASTQ files did not match the R1/R2 (or _1/_2) pattern: {_preview(unmatched)}")
    if duplicates:
        warn(f"{len(duplicates)} FASTQ files duplicate an already paired read and were ignored: {_preview(duplicates)}")
    if orphans:
        warn(f"{len(orphans)} samples missing R1 or R2: {_preview(orphans)}")
    return {s: p for s, p in pairs.items() if len(p) == 2}


def read_samplesheet(path):
    """Read a CSV/TSV with columns sample,r1,r2. Relative read paths are resolved against the sheet's folder."""
    path = os.path.abspath(os.path.expanduser(path))
    with open(path, newline="") as fh:
        lines = [line for line in fh.read().splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"samplesheet is empty: {path}")
    dialect = "excel-tab" if "\t" in lines[0] else "excel"
    pairs = {}
    for row in csv.DictReader(lines, dialect=dialect):
        row = {(k or "").strip().lower(): (v or "").strip() for k, v in row.items()}
        sample = row.get("sample", "")
        if not SAMPLE_NAME_RE.match(sample):
            raise ValueError(f"invalid sample name '{sample}' (use letters, digits, '.', '_' or '-')")
        if sample in pairs:
            raise ValueError(f"duplicate sample '{sample}' in samplesheet")
        pairs[sample] = {}
        for key in ("r1", "r2"):
            fastq = os.path.join(os.path.dirname(path), os.path.expanduser(row.get(key, "")))
            if not row.get(key) or not os.path.isfile(fastq):
                raise ValueError(f"{key.upper()} for sample '{sample}' not found: '{row.get(key, '')}'")
            pairs[sample][key.upper()] = fastq
    ids = [fastq_id(f) for reads in pairs.values() for f in reads.values()]
    clashes = sorted({i for i in ids if i and ids.count(i) > 1})
    if None in ids or clashes:
        raise ValueError(f"read files must end in {'/'.join(FASTQ_EXTENSIONS)} and have unique names: {clashes}")
    return pairs


def pangenome_groups(results_csv):
    """{species: [samples]} for every assigned species with at least two samples."""
    import pandas as pd

    df = pd.read_csv(results_csv, dtype=str, keep_default_na=False)
    groups = {}
    for species, samples in df.groupby("SPECIES")["SAMPLE"]:
        unique = sorted(set(samples))
        if species not in UNASSIGNED_SPECIES and len(unique) > 1:
            groups[species] = unique
    return groups


def mlst_scheme(mlst_tsv):
    """The MLST scheme name (column 2 of tseemann/mlst output), or '' if unavailable."""
    try:
        with open(mlst_tsv) as fh:
            fields = fh.readline().rstrip("\n").split("\t")
    except OSError:
        return ""
    return fields[1] if len(fields) > 1 else ""


def amrfinder_organism(mlst_tsv):
    """AMRFinderPlus --organism value for a sample, or None when the species has no point-mutation set."""
    name = species_dict.get(mlst_scheme(mlst_tsv), "")
    return AMRFINDER_ORGANISMS.get(" ".join(name.split()[:2]))
