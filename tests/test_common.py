import os

import pandas as pd
import pytest

from pegas import common
from pegas.common import (
    amrfinder_organism,
    build_fastq_pairs,
    fastq_id,
    pangenome_groups,
    parse_fastq_read,
    read_samplesheet,
)


@pytest.mark.parametrize("name, expected", [
    ("001_R1.fastq.gz", ("001", "R1")),
    ("S1_L001_R2_001.fastq.gz", ("S1_L001", "R2")),
    ("CRE25-RO32-EC26_S12_ME_R1_L001.fastq.gz", ("CRE25-RO32-EC26_S12_ME", "R1")),
    ("sampleR2.fq.gz", ("sample", "R2")),
    ("SRR12345_1.fastq.gz", ("SRR12345", "R1")),
    ("isolate-2_R1.fq", ("isolate-2", "R1")),
    ("strain_12.fastq.gz", None),
    ("SRR12_x.fastq.gz", None),
    ("reads_R1.bam", None),
])
def test_parse_fastq_read(name, expected):
    assert parse_fastq_read(name) == expected


def test_fastq_id_matches_fastqc_naming():
    assert fastq_id("/x/a_R1.fq.gz") == "a_R1"
    assert fastq_id("a_R1.fastq") == "a_R1"
    assert fastq_id("a_R1.txt") is None


def test_build_fastq_pairs_drops_orphans_and_unmatched(capsys):
    pairs = build_fastq_pairs(["/d/A_R1.fastq.gz", "/d/A_R2.fastq.gz", "/d/B_1.fq.gz", "/d/odd.fastq.gz"])
    assert pairs == {"A": {"R1": "/d/A_R1.fastq.gz", "R2": "/d/A_R2.fastq.gz"}}
    err = capsys.readouterr().err
    assert "odd.fastq.gz" in err and "missing R1 or R2: B" in err


def test_read_samplesheet_resolves_relative_paths(tmp_path):
    for name in ("x_fwd.fq.gz", "x_rev.fq.gz"):
        (tmp_path / name).write_text("")
    sheet = tmp_path / "samples.csv"
    sheet.write_text("sample,r1,r2\nX,x_fwd.fq.gz,x_rev.fq.gz\n")
    assert read_samplesheet(sheet) == {"X": {"R1": str(tmp_path / "x_fwd.fq.gz"), "R2": str(tmp_path / "x_rev.fq.gz")}}


@pytest.mark.parametrize("content, message", [
    ("sample,r1,r2\nX,missing.fq.gz,x.fq.gz\n", "not found"),
    ("sample,r1,r2\nbad name,x.fq.gz,x.fq.gz\n", "invalid sample name"),
    ("sample\tr1\tr2\nX\tx.fq.gz\tx.fq.gz\nX\tx.fq.gz\tx.fq.gz\n", "duplicate sample"),
])
def test_read_samplesheet_errors(tmp_path, content, message):
    (tmp_path / "x.fq.gz").write_text("")
    sheet = tmp_path / "s.tsv"
    sheet.write_text(content)
    with pytest.raises(ValueError, match=message):
        read_samplesheet(sheet)


@pytest.mark.parametrize("scheme, organism", [
    ("ecoli_achtman_4", "Escherichia"),
    ("klebsiella", "Klebsiella_pneumoniae"),
    ("abaumannii_2", "Acinetobacter_baumannii"),
    ("vcholerae_2", "Vibrio_cholerae"),
    ("neisseria", None),
    ("-", None),
])
def test_amrfinder_organism(tmp_path, scheme, organism):
    mlst = tmp_path / "mlst.tsv"
    mlst.write_text(f"contigs.fa\t{scheme}\t131\tadk(1)\n")
    assert amrfinder_organism(mlst) == organism


def test_amrfinder_organism_without_mlst(tmp_path):
    assert amrfinder_organism(tmp_path / "absent.tsv") is None


def test_pangenome_groups_skip_unassigned_and_singletons(tmp_path):
    csv = tmp_path / "results.csv"
    pd.DataFrame({
        "SAMPLE": ["a", "a", "b", "c", "d", "e", "f"],
        "SPECIES": ["ecoli", "ecoli", "ecoli", "saureus", "Unassigned", "Unassigned", "-"],
    }).to_csv(csv, index=False)
    assert pangenome_groups(csv) == {"ecoli": ["a", "b"]}


def test_clean_env_drops_foreign_libraries(monkeypatch):
    monkeypatch.setenv("PERL5LIB", "/home/user/perl5/lib")
    monkeypatch.setenv("PYTHONPATH", "/opt/other")
    env = common.clean_env()
    assert "PERL5LIB" not in env and "PYTHONPATH" not in env
    assert env["PYTHONNOUSERSITE"] == "1" and env["PATH"] == os.environ["PATH"]
