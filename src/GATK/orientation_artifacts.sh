#!/bin/bash

# GATK GetPileupSummaries script for orientation artifacts

# Input: BAM file
# Output: Summary of orientation artifacts

INPUT_BAM=$1  # Input BAM file
OUTPUT_FILE=$2  # Output summary file

# Run GetPileupSummaries
gatk GetPileupSummaries \
  -I $INPUT_BAM \
  -R $REFERENCE_GENOME \
  -O $OUTPUT_FILE

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "GetPileupSummaries failed"
  exit 1
else
  echo "GetPileupSummaries completed successfully"
fi
