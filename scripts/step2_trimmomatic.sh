#!/bin/bash

# Step 2: Trim Illumina paired-end reads using Trimmomatic

# Load Trimmomatic module
module load trimmomatic/0.39-tfj436w

# Process paired-end reads
for file in *_R1_cut.fastq
do
    withpath="${file}"
    filename=${withpath##*/}  
    base="${filename%_R1_cut.fastq}"
    echo "${base}"

    # Run Trimmomatic
    java -jar /path/to/trimmomatic.jar PE -phred33 -trimlog trimmomatic_summary.txt \
        "${base}_R1_cut.fastq" \
        "${base}_R2_cut.fastq" \
        "${base}_Q15_R1_cut.trimmed_paired.fastq" "${base}_Q15_R1_cut.trimmed_unpaired.fastq" \
        "${base}_Q15_R2_cut.trimmed_paired.fastq" "${base}_Q15_R2_cut.trimmed_unpaired.fastq" \
        ILLUMINACLIP:/path/to/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
