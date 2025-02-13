# conda install -c bioconda trimmomatic -y
# mkdir -p adapters
# wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa -O adapters/TruSeq3-PE.fa
# Create directory for trimmed files
mkdir trimmed_reads

# Trim Cont1
# Trimmomatic Wall time: 19min 37s, Wall time: 20min 45s
trimmomatic PE \
  Tien_Cont1_S25_R1_001.fastq.gz Tien_Cont1_S25_R2_001.fastq.gz \
  trimmed_reads/Cont1_R1_paired.fastq.gz trimmed_reads/Cont1_R1_unpaired.fastq.gz \
  trimmed_reads/Cont1_R2_paired.fastq.gz trimmed_reads/Cont1_R2_unpaired.fastq.gz \
  ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36

# Trim Zn2
trimmomatic PE \
  Tien_Zn2_S29_R1_001.fastq.gz Tien_Zn2_S29_R2_001.fastq.gz \
  trimmed_reads/Zn2_R1_paired.fastq.gz trimmed_reads/Zn2_R1_unpaired.fastq.gz \
  trimmed_reads/Zn2_R2_paired.fastq.gz trimmed_reads/Zn2_R2_unpaired.fastq.gz \
  ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36