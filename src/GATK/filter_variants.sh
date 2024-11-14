#!/bin/bash

# GATK VariantFiltration script

# Input: VCF file
# Output: Filtered VCF file

INPUT_VCF=$1  # Input VCF file
OUTPUT_VCF=$2  # Output filtered VCF file

# Run VariantFiltration
gatk VariantFiltration \
  -R $REFERENCE_GENOME \
  -V $INPUT_VCF \
  --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 60.0" \
  -O $OUTPUT_VCF

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "VariantFiltration failed"
  exit 1
else
  echo "VariantFiltration completed successfully"
fi
