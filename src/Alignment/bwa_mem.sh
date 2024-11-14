#!/bin/bash

# BWA MEM alignment script

# Input: Reference genome, paired-end reads
# Output: SAM file

REF_GENOME=$1  # Path to reference genome
R1=$2  # Read 1 file
R2=$3  # Read 2 file
OUTPUT=$4  # Output SAM file

# Run BWA MEM
bwa mem -t 4 $REF_GENOME $R1 $R2 > $OUTPUT

# Check if the alignment was successful
if [ $? -ne 0 ]; then
  echo "BWA MEM failed"
  exit 1
else
  echo "BWA MEM completed successfully"
fi
