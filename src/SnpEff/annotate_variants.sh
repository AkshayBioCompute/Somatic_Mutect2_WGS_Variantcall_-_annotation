#!/bin/bash

# SnpEff annotation script

# Input: VCF file
# Output: Annotated VCF file

INPUT_VCF=$1  # Input VCF file
OUTPUT_VCF=$2  # Output annotated VCF file

# Run SnpEff annotation
java -Xmx4g -jar snpEff.jar -v GRCh38.99 $INPUT_VCF > $OUTPUT_VCF

# Check if the command was successful
if [ $? -ne 0 ]; then
  echo "SnpEff annotation failed"
  exit 1
else
  echo "SnpEff annotation completed successfully"
fi
