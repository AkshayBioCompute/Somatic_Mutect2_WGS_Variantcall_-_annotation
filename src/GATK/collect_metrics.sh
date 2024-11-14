#!/bin/bash

# GATK CollectAlignmentSummaryMetrics script

# Input: BAM file
# Output: Alignment summary metrics

INPUT_BAM=$1  # Input BAM file
OUTPUT_FILE=$2  # Output metrics file

# Run CollectAlignmentSummaryMetrics
gatk CollectAlignmentSummaryMetrics \
  -I $INPUT_BAM \
  -R $REFERENCE_GENOME \
  -O $OUTPUT_FILE

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "CollectAlignmentSummaryMetrics failed"
  exit 1
else
  echo "CollectAlignmentSummaryMetrics completed successfully"
fi
