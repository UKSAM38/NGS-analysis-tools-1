#!/bin/bash

# Step 4: Analyse BAM files using GATK

# Activate GATK environment
conda activate gatk_env

# Create output directory
mkdir -p ASM_output

# Analyze BAM files using GATK
for file in *_sorted.bam
do
    withpath="${file}"
    filename=${withpath##*/}  
    base="${filename%_sorted.bam}"
    echo "Analyzing ${base}"

    # Run GATK AnalyzeSaturationMutagenesis
    /path/to/gatk/gatk AnalyzeSaturationMutagenesis -I "${file}" -R ReferenceGenome.fasta --orf 106-846 -O ASM_output/"${base}_Q15_output.bam"
done
