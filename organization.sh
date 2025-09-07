#!/bin/bash

# Usage: `bash organization.sh ./Figures png Outputs`
# Target directory for new files
target_dir="$1"
extension="$2"
source_dir="$3"

# Ensure the target directory exists
mkdir -p "$target_dir"

# Use find to locate files and copy them with transformation
find $source_dir -type f -name "*.$extension" | while read file; do
    # Transform the filename by replacing '/' with '-'
    new_filename=$(echo "$file" | sed 's|/|-|g')
    
    # Formulate the full path for the new file
    cp "$file" "$target_dir/$new_filename"
    
    # Print a confirmation or debug message (optional)
    echo "Copied $file to $target_dir/$new_filename"
done
