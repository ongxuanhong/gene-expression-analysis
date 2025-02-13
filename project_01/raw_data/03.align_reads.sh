# First, build index if not already available
# hisat2-build reference_genome.fa genome_index

# Create directory for aligned files
mkdir -p aligned_reads

for sample in Tien_Cont1_S25 Tien_Cont2_S26 Tien_Cont3_S27 Tien_Zn1_S28 Tien_Zn2_S29 Tien_Zn3_S30
do
    hisat2 -p 8 \
        -x hisat2_index/genome_index \
        -1 trimmed_reads/${sample}_R1_paired.fastq.gz \
        -2 trimmed_reads/${sample}_R2_paired.fastq.gz \
        -S aligned_reads/${sample}.sam
done

# Create directory for processed BAM files
mkdir -p bam_files

for sample in Tien_Cont1_S25 Tien_Cont2_S26 Tien_Cont3_S27 Tien_Zn1_S28 Tien_Zn2_S29 Tien_Zn3_S30
do
    # Convert SAM to BAM
    samtools view -bS aligned_reads/${sample}.sam > bam_files/${sample}.bam
    
    # Sort BAM file
    samtools sort bam_files/${sample}.bam -o bam_files/${sample}.sorted.bam
    
    # Index BAM file
    samtools index bam_files/${sample}.sorted.bam
    
    # Remove intermediate files
    rm aligned_reads/${sample}.sam bam_files/${sample}.bam
done