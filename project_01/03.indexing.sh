# conda install -c bioconda star -y
# download reference and annotation files
# Genome (DNA) FASTA
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.transcripts.fa.gz
# GTF annotation
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh37_mapping/gencode.v47lift37.annotation.gtf.gz
# For hg38 (GRCh38) - Release 44
# Genome sequence
mkdir reference && cd reference

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
# Gene annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

# For hg19 (GRCh37) - Release 19
# Genome sequence
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
# Gene annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

# Decompress all files
gunzip *.gz

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
# HISAT Wall time: 14min 14s
hisat2-build gencode.v47lift37.transcripts.fa hisat2_index/genome_index

# Worked
# HISAT alignment Wall time: 1h 3min 27s Wall time: 1h 2min 3s
hisat2 -p 8 \
-x hisat2_index/genome_index \
-1 trimmed_reads/Cont1_R1_paired.fastq.gz \
-2 trimmed_reads/Cont1_R2_paired.fastq.gz \
| samtools view -bS - \
| samtools sort - -o alignments/Cont1.sorted.bam  

hisat2 -p 8 \
-x hisat2_index/genome_index \
-1 trimmed_reads/Zn2_R1_paired.fastq.gz \
-2 trimmed_reads/Zn2_R2_paired.fastq.gz \
| samtools view -bS - \
| samtools sort - -o alignments/Zn2.sorted.bam  


# Index the BAM files
samtools index alignments/Cont1.sorted.bam
samtools index alignments/Zn2.sorted.bam