# Create directories to organize your analysis
mkdir -p fastqc_results
mkdir -p multiqc_results

# Run FastQC on all fastq.gz files
# conda install -c bioconda fastqc multiqc -y
fastqc \
  raw_data/Tien_Cont1_S25_R1_001.fastq.gz \
  --outdir fastqc_results \
  --threads 4

# Alternative: you can use wildcard to run all fastq.gz files at once
# fastqc *.fastq.gz --outdir fastqc_results --threads 4

# Run MultiQC to aggregate FastQC results
multiqc fastqc_results/ -o multiqc_results/