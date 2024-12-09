# Guide for NGS Data Analysis

This guide outlines a 4-step process for analyzing NGS data using the provided scripts. Each step includes instructions and necessary prerequisites.

---

## Step 1: Remove Primer Sequences with Cutadapt (PCR1)

1. **Set the Working Directory**  
   Ensure the working directory is correctly set in the script:
   ```bash
   cd /path/to/working/directory
   ```

2. **Activate Conda Environment**  
   Use the appropriate conda environment for `cutadapt`:
   ```bash
   eval "$(conda shell.bash hook)"
   conda activate cutadapt_env
   ```

3. **Run the Script**  
   Execute the `step1_cutadapt.sh` script to remove primer sequences from paired reads:
   ```bash
   bash scripts/step1_cutadapt.sh
   ```

4. **Output Files**  
   - Processed files are saved as `_R1_cut.fastq` and `_R2_cut.fastq`.

---

## Step 2: Trim Illumina Paired-End Reads with Trimmomatic

1. **Load Trimmomatic Module**  
   Ensure Trimmomatic is installed and loaded:
   ```bash
   module load trimmomatic/0.39-tfj436w
   ```

2. **Run the Script**  
   Execute the `step2_trimmomatic.sh` script to trim low-quality reads:
   ```bash
   bash scripts/step2_trimmomatic.sh
   ```

3. **Output Files**  
   - Trimmed paired and unpaired reads are saved with `_trimmed_paired.fastq` and `_trimmed_unpaired.fastq` suffixes.

---

## Step 3: Align Paired Reads Using BWA and Create BAM Files

1. **Prepare Reference Genome**  
   Ensure the reference genome is indexed:
   ```bash
   bwa index -a is ReferenceGenome.fasta
   samtools faidx ReferenceGenome.fasta
   java -jar /path/to/picard.jar CreateSequenceDictionary R=ReferenceGenome.fasta O=ReferenceGenome.dict
   ```

2. **Run the Script**  
   Execute the `step3_bwa_samtools.sh` script to align paired reads and create sorted BAM files:
   ```bash
   bash scripts/step3_bwa_samtools.sh
   ```

3. **Output Files**  
   - Sorted BAM files are saved with `_sorted.bam` suffixes.

---

## Step 4: Analyze BAM Files Using GATK

1. **Activate GATK Environment**  
   Use the appropriate conda environment for GATK:
   ```bash
   conda activate gatk_env
   ```

2. **Create Output Directory**  
   Ensure the directory for GATK output exists:
   ```bash
   mkdir -p ASM_output
   ```

3. **Run the Script**  
   Execute the `step4_gatk.sh` script to perform variant analysis:
   ```bash
   bash scripts/step4_gatk.sh
   ```

4. **Output Files**  
   - Results are saved in the `ASM_output` directory with analysis-specific suffixes.

---

## Notes

- Replace placeholder paths (e.g., `/path/to/working/directory`) with your actual directory structure.
- Ensure all required tools (Cutadapt, Trimmomatic, BWA, SAMtools, GATK) are installed and accessible.
- Adjust parameters in the scripts to suit your dataset and analysis goals.

