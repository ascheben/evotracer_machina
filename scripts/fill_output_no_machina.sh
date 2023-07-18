#!/bin/bash

# Input and output file paths
input_file=""
output_file=""

# Process the CSV file using Awk
awk -F ',' -v OFS=',' '{
  # Get the number of fields in the line
  num_fields = NF

  # Get the values of the 3rd to last, 2nd to last, and last fields
  third_last = $(num_fields - 2)
  second_last = $(num_fields - 1)
  last = $num_fields

  # Remove leading/trailing spaces from the fields
  gsub(/^[[:blank:]]+|[[:blank:]]+$/, "", third_last)
  gsub(/^[[:blank:]]+|[[:blank:]]+$/, "", second_last)
  gsub(/^[[:blank:]]+|[[:blank:]]+$/, "", last)

  # Check the specified conditions and update the values if necessary
  if (third_last == "0") {
    second_last = "0"
    last = "1.0000"
  }

if (third_last != "0" && last == "") {
    second_last = "0"
    last = "0.0000"
  }

  # Print the modified line or the original line if no modification is made
  for (i = 1; i <= num_fields - 3; i++) {
    printf "%s%s", $i, OFS
  }
  print third_last, second_last, last
}' "$input_file" > "$output_file"
