#!/bin/bash

# Usage:
#   ./cut_all_stl.sh INPUT_DIR OUTPUT_DIR reference.stl
#
# Example:
#   ./cut_all_stl.sh ./meshes ./cut_meshes ./reference.stl

set -u

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 INPUT_DIR OUTPUT_DIR REFERENCE_STL"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
REFERENCE_STL="$3"

# Location of intersect_stl.py.
# Change this if the Python script is somewhere else.
INTERSECT_SCRIPT="./intersect_stl.py"

# Check input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Check reference file
if [ ! -f "$REFERENCE_STL" ]; then
    echo "Error: reference STL does not exist: $REFERENCE_STL"
    exit 1
fi

# Check intersection script
if [ ! -f "$INTERSECT_SCRIPT" ]; then
    echo "Error: cannot find: $INTERSECT_SCRIPT"
    exit 1
fi

# Create output directory if necessary
mkdir -p "$OUTPUT_DIR"

count=0
success=0
failed=0

# Process all .stl files
for input_file in "$INPUT_DIR"/*.stl; do

    # If there are no matching files, stop
    [ -e "$input_file" ] || continue

    filename=$(basename "$input_file")
    basename="${filename%.*}"

    output_file="$OUTPUT_DIR/${basename}.stl"

    echo "------------------------------------------------------------"
    echo "Input:  $input_file"
    echo "Output: $output_file"

    count=$((count + 1))

    if python "$INTERSECT_SCRIPT" \
        "$input_file" \
        "$REFERENCE_STL" \
        "$output_file"
    then
        success=$((success + 1))
    else
        echo "WARNING: intersection failed or was empty for $filename"
        failed=$((failed + 1))
    fi

done

echo "------------------------------------------------------------"
echo "Finished."
echo "Files processed: $count"
echo "Successful:      $success"
echo "Failed/empty:    $failed"
