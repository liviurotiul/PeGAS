#!/bin/bash

echo "Starting custom build script..."

# Ensure the script is executable
chmod +x install.sh

# Define the destination directory
DEST_DIR="$PREFIX/share/pegas"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Copy the Snakefile and assets to the destination directory
cp src/pegas/Snakefile "$DEST_DIR"
cp src/pegas/layout.html "$DEST_DIR"
cp src/pegas/gc_content.json "$DEST_DIR"

# Copy conda envs into a dedicated folder
mkdir -p "$DEST_DIR/envs"
cp src/pegas/envs/*.yml "$DEST_DIR/envs/"
