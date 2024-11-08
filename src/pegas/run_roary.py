import os
import subprocess
import pandas as pd
from tqdm import tqdm

tqdm.write("[pegas]Running pangenome analysis")

# Create pangenome directory if it does not exist
os.makedirs('pangenome', exist_ok=True)

# Read eligible species from dataframe
df = pd.read_csv('dataframe/results.csv', dtype={'SAMPLE': str})
df_grouped = df.groupby('SPECIES')['SAMPLE'].nunique().reset_index()
eligible_species = df_grouped.loc[df_grouped['SAMPLE'] > 1, 'SPECIES'].tolist()

tqdm.write(f"[pegas]Eligible species: {eligible_species}")

# Iterate over each species and process
for species in eligible_species:
    tqdm.write(f"[pegas]Processing species: {species}")

    species_dir = os.path.join('pangenome', species)
    os.makedirs(species_dir, exist_ok=True)

    # Remove any existing output folder for species
    output_dir = os.path.join(species_dir, 'output')
    if os.path.exists(output_dir):
        subprocess.run(['rm', '-rf', output_dir])

    # Get associated samples for the species and copy GFF files
    samples = df.loc[df['SPECIES'] == species, 'SAMPLE'].unique().tolist()
    tqdm.write(f"[pegas]Samples for species {species}: {samples}")

    for sample in samples:

        src_gff = f"results/{sample}/prokka/{sample}.gff"
        dest_gff = os.path.join(species_dir, f"{sample}.gff")

        if not os.path.isfile(dest_gff):
            tqdm.write(f"[pegas]Copying {sample}.gff to pangenome/{species}")
            subprocess.run(['cp', src_gff, dest_gff])

    # Run roary for the species
    tqdm.write(f"[pegas]Running roary for species {species}")
    gff_files = [os.path.join(species_dir, f) for f in os.listdir(species_dir) if f.endswith('.gff')]

    if len(gff_files) > 1:
        roary_cmd = [
            'roary',
            '-f', output_dir,
            '-e',
            '-n',
            '-p', '32'
        ] + gff_files

        subprocess.run(roary_cmd)

        # Check if the output directory needs to be renamed
        for item in os.listdir(species_dir):
            if item.startswith('output_') and os.path.isdir(os.path.join(species_dir, item)):
                old_output = os.path.join(species_dir, item)
                subprocess.run(['mv', old_output, output_dir])
                tqdm.write(f"[pegas]Renaming {item} to output")
    else:
        tqdm.write(f"[pegas]Not enough samples for species {species} to run roary.")

tqdm.write("[pegas]Pangenome analysis finished, planting flag")

# Create the flag file
if not os.path.exists('flags'):
    os.makedirs('flags')

with open('flags/.pangenome', 'w') as f:
    f.write('')
