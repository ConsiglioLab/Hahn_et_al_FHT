#!/bin/bash
# This script copies all the salmon quant.sf files from a source directory to a destination directory.
# It renames the files by appending the parent sample name to the file name.

# Source directory where you want to search for files
source_directory="../runs/20240311_rnaseq39/salmon"

# Destination directory where you want to copy the files
destination_directory="../data/salmon_quant_files"

# if the destination directory does not exist, create it
if [ ! -d "$destination_directory" ]; then
    mkdir -p "$destination_directory"
fi

# Use the find command to search for specific files (e.g., all .txt files)
# Modify the find command with your desired criteria
find "$source_directory" -type f -name "quant.sf" -print0 |
while IFS= read -r -d '' file; do
    # Get the immediate parent directory name
    parent_dir=$(basename "$(dirname "$file")")

    # Extract the filename
    filename=$(basename "$file")

    # Copy the found files to the destination directory
    cp "$file" "$destination_directory/$parent_dir-$filename"
done

echo "Files copied successfully!"
