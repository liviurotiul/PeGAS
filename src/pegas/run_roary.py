"""Pangenome step for the Snakemake checkpoint: one Roary run per species group.

Usage: run_roary.py '<{"species": ["sample", ...]}>' <threads>

Runs as a plain script inside the Roary env (Python 3.7), so it must not use the
`snakemake` object. The groups come from pegas.common.pangenome_groups, so
samples without an assigned species are never pooled.
"""
import json
import os
import shutil
import subprocess
import sys

groups = json.loads(sys.argv[1])
threads = sys.argv[2]
os.makedirs("pangenome", exist_ok=True)
print("[pegas] Pangenome groups: {}".format(groups), flush=True)

for species, samples in groups.items():
    species_dir = os.path.join("pangenome", species)
    output_dir = os.path.join(species_dir, "output")
    os.makedirs(species_dir, exist_ok=True)
    shutil.rmtree(output_dir, ignore_errors=True)

    for sample in samples:
        src_gff = "results/{0}/prokka/{0}.gff".format(sample)
        dest_gff = os.path.join(species_dir, sample + ".gff")
        if not os.path.isfile(src_gff) or os.path.getsize(src_gff) == 0:
            print("[pegas] Skipping {} for Roary: valid Prokka GFF not found.".format(sample), flush=True)
            continue
        if os.path.lexists(dest_gff):
            os.remove(dest_gff)
        os.symlink(os.path.relpath(src_gff, species_dir), dest_gff)

    gff_files = sorted(os.path.join(species_dir, f) for f in os.listdir(species_dir) if f.endswith(".gff"))
    if len(gff_files) < 2:
        print("[pegas] Not enough samples for {} to run Roary.".format(species), flush=True)
        continue

    print("[pegas] Running Roary for {} ({} samples)".format(species, len(gff_files)), flush=True)
    subprocess.run(["roary", "-f", output_dir, "-e", "-n", "-p", threads] + gff_files)
    # Roary appends a timestamp when the output folder already exists.
    for item in os.listdir(species_dir):
        if item.startswith("output_") and os.path.isdir(os.path.join(species_dir, item)):
            shutil.move(os.path.join(species_dir, item), output_dir)
