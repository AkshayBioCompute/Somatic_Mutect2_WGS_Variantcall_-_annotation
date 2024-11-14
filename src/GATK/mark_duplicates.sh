#!/bin/bash

# GATK MarkDuplicatesSpark script

# Input: BAM file
# Output: BAM file with duplicates marked

INPUT_BAM=$1  # Input BAM file
OUTPUT_BAM=$2  # Output BAM file
METRICS_FILE=$3  # Output metrics file

# Run MarkDuplicatesSpark
gatk MarkDuplicatesSpark \
  -I $INPUT_BAM \
  -O $OUTPUT_BAM \
  --CREATE_INDEX true \
  --METRICS_FILE $METRICS_FILE

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "MarkDuplicatesSpark failed"
  exit 1
else
  echo "MarkDuplicatesSpark completed successfully"
fi
