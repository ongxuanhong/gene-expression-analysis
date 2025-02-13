# Create directory for trimmed files
mkdir trimmed_reads

for sample in Tien_Cont1_S25 Tien_Cont2_S26 Tien_Cont3_S27 Tien_Zn1_S28 Tien_Zn2_S29 Tien_Zn3_S30
do
    trimmomatic PE \
        ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz \
        trimmed_reads/${sample}_R1_paired.fastq.gz trimmed_reads/${sample}_R1_unpaired.fastq.gz \
        trimmed_reads/${sample}_R2_paired.fastq.gz trimmed_reads/${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
done