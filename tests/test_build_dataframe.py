import json
import os
import shutil
from pathlib import Path

import pandas as pd
import pytest

from pegas.build_dataframe import assess_qc, build_dataframe
from pegas.common import pangenome_groups

DATA = Path(__file__).parent / "data"


def write_sample(root, sample, scheme, st, amrfinder=None, plasmidfinder=None):
    sample_dir = root / "results" / sample
    (sample_dir / "shovill").mkdir(parents=True)
    contig = "ACGTACGT" * 2000
    (sample_dir / "shovill" / "contigs.fa").write_text(f">contig00001 len={len(contig)} cov=30.0\n{contig}\n")
    (sample_dir / "mlst.tsv").write_text(f"contigs.fa\t{scheme}\t{st}\n")
    for source, name in ((amrfinder, "amrfinder.tsv"), (plasmidfinder, "abricate_plasmidfinder.tsv")):
        if source:
            shutil.copy(source, sample_dir / name)
        else:
            (sample_dir / name).write_text("")


@pytest.fixture
def run_dir(tmp_path, monkeypatch):
    write_sample(tmp_path, "EC1", "ecoli_achtman_4", "131", DATA / "amrfinder_ecoli.tsv", DATA / "plasmidfinder_ecoli.tsv")
    write_sample(tmp_path, "EC2", "ecoli_achtman_4", "10", DATA / "amrfinder_ecoli.tsv")
    write_sample(tmp_path, "X1", "-", "-")
    pairs = {s: {"R1": f"/in/{s}_R1.fastq.gz", "R2": f"/in/{s}_R2.fastq.gz"} for s in ("EC1", "EC2", "X1")}
    (tmp_path / "sample_manifest.json").write_text(json.dumps({"sample_pairs": pairs}))
    monkeypatch.chdir(tmp_path)
    build_dataframe()
    return tmp_path


def test_amrfinder_hits_are_mapped(run_dir):
    df = pd.read_csv("dataframe/results.csv", keep_default_na=False)
    ec1 = df[df.SAMPLE == "EC1"]
    assert set(ec1.PREDICTION_SOURCE) == {"AMR", "Virulence", "Stress", "PlasmidFinder"}
    gyra = ec1[ec1.GENE == "gyrA_S83L"].iloc[0]
    assert (gyra.ELEMENT_SUBTYPE, gyra.RESISTANCE, gyra.SEQUENCE) == ("POINT", "QUINOLONE", "contig00002")
    assert set(ec1[ec1.PREDICTION_SOURCE == "Virulence"].RESISTANCE) == {""}


def test_unassigned_species_is_not_pooled(run_dir):
    df = pd.read_csv("dataframe/results.csv", keep_default_na=False)
    assert df[df.SAMPLE == "X1"].SPECIES.tolist() == ["Unassigned"]
    assert pangenome_groups("dataframe/results.csv") == {"ecoli_achtman_4": ["EC1", "EC2"]}


def test_status_and_qc_tables(run_dir):
    status = pd.read_csv("dataframe/sample_status.tsv", sep="\t", index_col="SAMPLE")
    assert status.loc["EC1", ["AMRFINDER", "PLASMIDFINDER", "MLST"]].tolist() == ["ok", "ok", "ok"]
    assert status.loc["X1", ["AMRFINDER", "MLST"]].tolist() == ["no_hits", "no_call"]
    qc = pd.read_csv("dataframe/assembly_qc.tsv", sep="\t", index_col="SAMPLE")
    assert str(qc.loc["EC1", "ST"]) == "131"
    # One 16 kb toy contig: N50 below 20 kb warns; 50% GC is fine for E. coli.
    assert (qc.loc["EC1", "QC_STATUS"], qc.loc["EC1", "QC_NOTES"]) == ("WARN", "n50 16000 < 20000")
    assert "species unassigned" in qc.loc["X1", "QC_NOTES"]


def test_presence_absence_matrices(run_dir):
    amr = pd.read_csv("dataframe/amr_matrix.tsv", sep="\t", index_col="SAMPLE")
    assert amr.loc["EC1", "gyrA_S83L"] == 1 and amr.loc["X1"].sum() == 0
    plasmid = pd.read_csv("dataframe/plasmid_matrix.tsv", sep="\t", index_col="SAMPLE")
    assert plasmid.loc["EC1"].sum() == 2 and plasmid.loc["EC2"].sum() == 0
    if os.path.exists("dataframe/pegas_results.xlsx"):
        sheets = pd.read_excel("dataframe/pegas_results.xlsx", sheet_name=None)
        assert {"sample_status", "assembly_qc", "amr_matrix"} <= set(sheets)


@pytest.mark.parametrize("row, verdict, note", [
    ({"SPECIES": "ecoli", "TOTAL_LENGTH": 5e6, "N50": 150000, "CONTIGS": 90, "GC_PERCENT": 50.6,
      "COMPLETENESS": 99.9, "CONTAMINATION": 0.4}, "PASS", ""),
    ({"SPECIES": "ecoli", "TOTAL_LENGTH": 5e6, "N50": 150000, "COMPLETENESS": 85}, "WARN", "completeness 85 < 90"),
    ({"SPECIES": "ecoli", "TOTAL_LENGTH": 5e6, "CONTAMINATION": 12}, "FAIL", "contamination 12 > 10"),
    ({"SPECIES": "ecoli", "TOTAL_LENGTH": 5e6, "GC_PERCENT": 39.0}, "WARN", "GC 39% vs expected 50.8%"),
    ({"SPECIES": "ecoli", "TOTAL_LENGTH": None}, "FAIL", "no assembly"),
])
def test_assess_qc(row, verdict, note):
    assert assess_qc(row) == (verdict, note)
