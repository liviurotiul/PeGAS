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
cp layout.html "$DEST_DIR"
