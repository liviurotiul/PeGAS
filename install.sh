#!/bin/bash

echo "Starting custom build script..."

# Ensure the script is executable
chmod +x install.sh

# Define the destination directory
DEST_DIR="$PREFIX/share/pegas"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Copy the Snakefile to the destination directory
cp src/pegas/Snakefile "$DEST_DIR"
cp src/pegas/layout.html "$DEST_DIR"
cp src/pegas/prokka_env.yml "$DEST_DIR"
cp src/pegas/mlst_env.yml "$DEST_DIR"
cp src/pegas/abricate_env.yml "$DEST_DIR"
cp src/pegas/shovill_env.yml "$DEST_DIR"
cp src/pegas/fastqc_env.yml "$DEST_DIR"
