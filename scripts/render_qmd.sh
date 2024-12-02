#!/bin/bash

# Read the qmd template file name & output document type from the command line
qmd_file=$1
doc_type=$2

# Get the qmd file name without the extension
qmd_name=${qmd_file%.*}

# Use jinja to create an intermediate qmd file
python ../../scripts/jinja2_apply.py $qmd_file $doc_type
# jinja $qmd_file -D output_format $doc_type -o $doc_type.qmd

# Render the intermediate qmd file to the desired format
quarto render $doc_type.qmd --to $doc_type -o ${qmd_name}.${doc_type}

# Remove the intermediate qmd file
rm $doc_type.qmd