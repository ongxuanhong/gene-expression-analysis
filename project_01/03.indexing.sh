# conda install -c bioconda star -y
# download reference and annotation files
# Genome (DNA) FASTA
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.transcripts.fa.gz
# GTF annotation
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.annotation.gtf.gz

# Uncompress
gunzip gencode.v47lift37.transcripts.fa.gz &&
gunzip gencode.v47lift37.annotation.gtf.gz

# https://github.com/alexdobin/STAR
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz &&
tar -xzf 2.7.11b.tar.gz &&
cd STAR-2.7.11b

# First create genome index
STAR --runMode genomeGenerate \
--genomeDir reference_genome_index \
--genomeFastaFiles gencode.v47lift37.transcripts.fa \
--sjdbGTFfile gencode.v47lift37.annotation.gtf \
--sjdbOverhang 99 \
--genomeSAindexNbases 12 \
--limitGenomeGenerateRAM 30000000000
# BLOCKED: EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM=30000000000is too small for your genome
# SOLUTION: please specify --limitGenomeGenerateRAM not less than 266994712842 and make that much RAM available 


# Create directory for index and alignments
mkdir hisat2_index
mkdir alignments

# Build HISAT2 index
# conda install hisat2 -y
hisat2-build gencode.v47lift37.transcripts.fa hisat2_index/genome_index

# Align Cont1
hisat2 -p 8 \
  -x hisat2_index/genome_index \
  -1 trimmed_reads/Cont1_R1_paired.fastq.gz \
  -2 trimmed_reads/Cont1_R2_paired.fastq.gz \
  -S alignments/Cont1.sam \
  --summary-file alignments/Cont1_alignment_summary.txt

# Align Zn2
hisat2 -p 8 \
  -x hisat2_index/genome_index \
  -1 trimmed_reads/Zn2_R1_paired.fastq.gz \
  -2 trimmed_reads/Zn2_R2_paired.fastq.gz \
  -S alignments/Zn2.sam \
  --summary-file alignments/Zn2_alignment_summary.txt

# Convert SAM to BAM, sort, and index
for sample in Cont1 Zn2
do
    samtools view -bS alignments/${sample}.sam > alignments/${sample}.bam
    samtools sort alignments/${sample}.bam -o alignments/${sample}.sorted.bam
    samtools index alignments/${sample}.sorted.bam
    rm alignments/${sample}.sam  # Remove SAM file to save space
done