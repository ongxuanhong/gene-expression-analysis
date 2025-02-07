# conda install -c bioconda trimmomatic -y
# For Cont1
trimmomatic PE \
Tien_Cont1_S25_R1_001.fastq.gz Tien_Cont1_S25_R2_001.fastq.gz \
Cont1_R1_paired.fastq.gz Cont1_R1_unpaired.fastq.gz \
Cont1_R2_paired.fastq.gz Cont1_R2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# For Zn2
trimmomatic PE \
Tien_Zn2_S29_R1_001.fastq.gz Tien_Zn2_S29_R2_001.fastq.gz \
Zn2_R1_paired.fastq.gz Zn2_R1_unpaired.fastq.gz \
Zn2_R2_paired.fastq.gz Zn2_R2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36