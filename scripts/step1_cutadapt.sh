#!/bin/bash

# Step 1: Remove primer sequences with cutadapt (PCR1)

# Navigate to the working directory
cd /path/to/working/directory

# Activate conda environment for cutadapt
eval "$(conda shell.bash hook)"
conda activate cutadapt_env

# Log output
exec > >(tee -a PCR1_cut.log) 2>&1

# Process each sample
for file in /path/to/working/directory/sample_R1.fastq.gz
do
    withpath="${file}"
    filename=${withpath##*/}  
    base="${filename%sample_R1.fastq.gz}"
    echo "${base}"

    # Run cutadapt
    cutadapt -e 0.2 -g ForwardPrimerSequence -G ReversePrimerSequence \
        -o "/path/to/working/directory/${base}_R1_cut.fastq" \
        -p "/path/to/working/directory/${base}_R2_cut.fastq" \
        "/path/to/working/directory/${base}sample_R1.fastq.gz" \
        "/path/to/working/directory/${base}sample_R2.fastq.gz"
done
