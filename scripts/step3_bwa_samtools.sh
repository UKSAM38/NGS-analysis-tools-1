#!/bin/bash

# Step 3: Align paired reads using BWA and create BAM files

# Load necessary HPC modules
module load bwa
module load samtools
module load picard

# Index reference genome
bwa index -a is ReferenceGenome.fasta
samtools faidx ReferenceGenome.fasta
java -jar /path/to/picard.jar CreateSequenceDictionary R=ReferenceGenome.fasta O=ReferenceGenome.dict

# Align paired reads and create BAM files
for file in *_R1_cut.trimmed_paired.fastq
do
    withpath="${file}"
    filename=${withpath##*/}  
    base="${filename%_R1_cut.trimmed_paired.fastq}"
    echo "Processing ${base}"

    # Align reads and generate SAM file
    bwa mem ReferenceGenome.fasta "${base}_R1_cut.trimmed_paired.fastq" "${base}_R2_cut.trimmed_paired.fastq" > "${base}.sam"

    # Convert SAM to BAM, sort and index
    samtools view -S -b "${base}.sam" -o "${base}.bam"
    samtools sort "${base}.bam" -o "${base}_sorted.bam"
    samtools index "${base}_sorted.bam"

    # Clean up
    rm "${base}.sam"
done
