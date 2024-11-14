#!/bin/bash

# GATK BaseRecalibrator script

# Input: BAM file, known sites for recalibration
# Output: Recalibration table

INPUT_BAM=$1  # Input BAM file
KNOWN_SITES=$2  # Known sites VCF file
OUTPUT_RECAL_TABLE=$3  # Output recalibration table

# Run BaseRecalibrator
gatk BaseRecalibrator \
  -I $INPUT_BAM \
  -R $REFERENCE_GENOME \
  --known-sites $KNOWN_SITES \
  -O $OUTPUT_RECAL_TABLE

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "BaseRecalibrator failed"
  exit 1
else
  echo "BaseRecalibrator completed successfully"
fi
