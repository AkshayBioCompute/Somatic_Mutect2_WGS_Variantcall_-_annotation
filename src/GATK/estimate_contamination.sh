#!/bin/bash

# GATK CalculateContamination script

# Input: VCF file
# Output: Contamination estimate

INPUT_VCF=$1  # Input VCF file
OUTPUT_FILE=$2  # Output contamination estimate file

# Run CalculateContamination
gatk CalculateContamination \
  -I $INPUT_VCF \
  -O $OUTPUT_FILE

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "CalculateContamination failed"
  exit 1
else
  echo "CalculateContamination completed successfully"
fi
