#!/bin/bash

# GATK ApplyBQSR script

# Input: BAM file, recalibration table
# Output: BAM file with applied BQSR

INPUT_BAM=$1  # Input BAM file
RECAL_TABLE=$2  # Recalibration table
OUTPUT_BAM=$3  # Output BAM file

# Run ApplyBQSR
gatk ApplyBQSR \
  -I $INPUT_BAM \
  -R $REFERENCE_GENOME \
  --bqsr $RECAL_TABLE \
  -O $OUTPUT_BAM

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "ApplyBQSR failed"
  exit 1
else
  echo "ApplyBQSR completed successfully"
fi
