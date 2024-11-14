#!/bin/bash

# BWA Index creation script

# Input: Reference genome
# Output: BWA index files (.amb, .ann, .bwt, .pac, .sa)

REF_GENOME=$1  # Path to reference genome

# Create the BWA index
bwa index $REF_GENOME

# Check if indexing was successful
if [ $? -ne 0 ]; then
  echo "BWA Index creation failed"
  exit 1
else
  echo "BWA Index creation completed successfully"
fi
