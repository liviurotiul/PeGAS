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
            (folder / f"{sample}_{mate}.fastq.gz").write_text("")
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
    for rule in ("shovill", "mlst", "amrfinder_organism", "amrfinder", "plasmidfinder", "build_dataframe"):
        assert rule in result.stdout
