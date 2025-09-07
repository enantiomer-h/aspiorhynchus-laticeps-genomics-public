#!/bin/bash

file_path="$1" # Replace with your file path
search_string="$2"

count=$(grep -o "$search_string" "$file_path" | wc -l)

echo "The string '$search_string' appears $count times in the file '$file_path'."
