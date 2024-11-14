#!/bin/bash

# GATK Mutect2 variant calling script

# Input: BAM files (tumor and normal)
# Output: VCF file with variants

TUMOR_BAM=$1  # Tumor BAM file
NORMAL_BAM=$2  # Normal BAM file
OUTPUT_VCF=$3  # Output VCF file

# Run Mutect2
gatk Mutect2 \
  -I $TUMOR_BAM \
  -I $NORMAL_BAM \
  -R $REFERENCE_GENOME \
  -O $OUTPUT_VCF

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "Mutect2 variant calling failed"
  exit 1
else
  echo "Mutect2 variant calling completed successfully"
fi
