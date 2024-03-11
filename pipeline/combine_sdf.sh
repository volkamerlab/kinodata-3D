#!/bin/bash

# Check if at least two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_folder> <output_file.sdf>"
    exit 1
fi

input_folder="$1"
output_file="$2"

# Check if input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Input folder does not exist: $input_folder"
    exit 1
fi

# # Check if output file already exists
# if [ -e "$output_file" ]; then
#     echo "Output file already exists: $output_file"
#     exit 1
# fi

# Initialize a temporary file for storing combined data
temp_file=$(mktemp)

# Loop through each input file in the folder
for input_file in "$input_folder"/*.sdf; do
    # Get filename without extension
    filename=$(basename -- "$input_file")
    filename="${filename%.*}"
    echo $filename

    # Add activity_id property to each entry in the SDF file
    head -n -1 $input_file >> $temp_file
    echo "> <activities.activity_id>" >> $temp_file
    echo "$filename" >> $temp_file
    echo "" >> $temp_file
    echo "\$\$\$\$" >> $temp_file
done

# Combine all data into a single SDF file
cat "$temp_file" > "$output_file"

# Clean up temporary file
rm "$temp_file"

echo "Combined SDF files in $input_folder into $output_file"

