# Find missing contigs.fasta files in the data directory and print the missing folder names.

import os

from tqdm import tqdm

# Path to the data directory.
data_dir = "results"

# List of all the folders in the data directory.
folders = os.listdir(data_dir)

# Loop through each folder.
for folder in folders:
    # Path to the contigs.fasta file.
    contigs_path = os.path.join(data_dir, folder, "shovill", "contigs.fasta")
    # If the file does not exist, print the name of the folder.
    if not os.path.exists(contigs_path):
        tqdm.write(folder)
