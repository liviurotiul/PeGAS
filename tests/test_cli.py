import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from pegas import main as cli

PEGAS_DIR = Path(cli.__file__).parent


def make_reads(folder, samples=("A", "B")):
    folder.mkdir(exist_ok=True)
    for sample in samples:
        for mate in ("R1", "R2"):
            with gzip.open(folder / f"{sample}_{mate}.fastq.gz", "wt") as fh:
                fh.write("@read1\nACGT\n+\nIIII\n")
    return folder


def test_prepare_output_dir_refuses_foreign_folder(tmp_path):
    (tmp_path / "thesis.docx").write_text("precious")
    with pytest.raises(SystemExit):
        cli.prepare_output_dir(str(tmp_path), overwrite=True)
    assert (tmp_path / "thesis.docx").exists()


def test_prepare_output_dir_modes(tmp_path):
    out = tmp_path / "out"
    assert cli.prepare_output_dir(str(out), overwrite=False) == "new"
    (out / cli.CONFIG_FILENAME).write_text("{}")
    (out / "results").mkdir()
    (out / ".snakemake").mkdir()
    assert cli.prepare_output_dir(str(out), overwrite=False) == "resume"
    assert (out / "results").exists()
    assert cli.prepare_output_dir(str(out), overwrite=True, dry_run=True) == "overwrite"
    assert (out / "results").exists()
    assert cli.prepare_output_dir(str(out), overwrite=True) == "overwrite"
    assert sorted(os.listdir(out)) == [".snakemake"]


def test_finalize_arguments_sets_every_core_option():
    for lite in (False, True):
        args = cli.finalize_arguments(cli.build_parser().parse_args(["-d", "x", "-o", "y", "-c", "8"]), lite_mode=lite)
        assert all(getattr(args, option) >= 1 for option in cli.CORE_OPTIONS)


def test_dry_run_lists_samples_without_touching_output(tmp_path, monkeypatch, capsys):
    reads = make_reads(tmp_path / "reads")
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["pegas", "-d", str(reads), "-o", str(out), "--dry-run"])
    cli.main()
    printed = capsys.readouterr().out
    assert "2 sample(s)" in printed and "A: A_R1.fastq.gz + A_R2.fastq.gz" in printed
    assert not out.exists()


def test_samples_with_empty_reads_are_skipped(tmp_path, monkeypatch, capsys):
    reads = make_reads(tmp_path / "reads", samples=("A", "B", "C"))
    with gzip.open(reads / "B_R2.fastq.gz", "wt"):
        pass
    monkeypatch.setattr(sys, "argv", ["pegas", "-d", str(reads), "-o", str(tmp_path / "out"), "--dry-run"])
    cli.main()
    captured = capsys.readouterr()
    assert "2 sample(s)" in captured.out and "B:" not in captured.out
    assert "Skipping B: no reads in B_R2.fastq.gz" in captured.err


def test_missing_input_is_reported(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["pegas", "-o", str(tmp_path / "out")])
    with pytest.raises(SystemExit):
        cli.main()


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_snakefile_dry_run(tmp_path):
    reads = make_reads(tmp_path / "reads")
    result = subprocess.run(
        ["snakemake", "-n", "--snakefile", str(PEGAS_DIR / "Snakefile"), "--directory", str(tmp_path / "out"),
         "--cores", "4", "--config", f"raw_data={reads}", f"install_path={PEGAS_DIR}",
         f"output_dir={tmp_path / 'out'}"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for rule in ("shovill", "mlst", "amrfinder_organism", "amrfinder", "abricate", "abricate_vfdb.tsv", "build_dataframe"):
        assert rule in result.stdout


def test_resume_keeps_samplesheet_samples(tmp_path):
    # Samplesheet names need not match file names: results/isolate_A must survive a resume.
    out = tmp_path / "out"
    for folder in ("isolate_A", "gone"):
        (out / "results" / folder).mkdir(parents=True)
    reads = ["/reads/A_forward.fq.gz", "/reads/A_reverse.fq.gz"]
    pairs = {"isolate_A": {"R1": reads[0], "R2": reads[1]}}
    cli.remove_extra_files(str(out), "raw_data", reads, pairs)
    assert sorted(os.listdir(out / "results")) == ["isolate_A"]


def test_lite_fastqc_failure_is_marked_on_the_samplesheet_sample(tmp_path, monkeypatch):
    from pegas import pipeline

    def broken(*args, **kwargs):
        raise subprocess.CalledProcessError(1, "fastqc")

    monkeypatch.setattr(pipeline, "run_command", broken)
    pairs = {"isolate_A": {"R1": "/reads/A_forward.fq.gz", "R2": "/reads/A_reverse.fq.gz"}}
    pipeline.run_fastqc(pairs, str(tmp_path), failure_policy="warn")
    assert (tmp_path / "results" / "isolate_A" / "status" / "fastqc.failed.txt").is_file()


def test_checkm2_db_dir_comes_from_the_environment(monkeypatch, tmp_path):
    import importlib
    from pegas import pipeline

    monkeypatch.setenv("CHECKM2_DB_DIR", str(tmp_path / "db"))
    try:
        assert importlib.reload(pipeline).CHECKM2_DB_DIR == str(tmp_path / "db")
    finally:
        monkeypatch.delenv("CHECKM2_DB_DIR")
        importlib.reload(pipeline)


def test_snakefile_fastqc_failure_is_handled():
    # Under --on-tool-failure warn a FastQC failure must mark the sample, not stop the run.
    text = (PEGAS_DIR / "Snakefile").read_text()
    rule = text.split("rule fastqc:")[1].split("\nrule ")[0]
    assert "fail {params.sample} fastqc {output}" in rule
