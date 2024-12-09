# Antibody Analysis Tools

This repository contains scripts for processing sequencing data to identify escape mutations in the receptor-binding domain (RBD).

## Workflow

1. **Step 1: Cutadapt**
   - Removes primer sequences from paired-end reads.

2. **Step 2: Trimmomatic**
   - Trims low-quality regions of the reads.

3. **Step 3: BWA and SAMTools**
   - Aligns paired-end reads and generates sorted BAM files.

4. **Step 4: GATK Analysis**
   - Identifies escape mutations using variant analysis.

## Requirements

- Cutadapt
- Trimmomatic
- BWA
- SAMTools
- GATK
- Python 3.8+

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/antibody-analysis-tools.git
   cd antibody-analysis-tools

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
