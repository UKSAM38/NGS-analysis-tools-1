#!/bin/bash
conda create -n ngs_analysis_env -y
conda activate ngs_analysis_env
conda install -c bioconda cutadapt trimmomatic bwa samtools gatk -y
echo "Environment setup complete!"
