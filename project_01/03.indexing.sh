# conda install -c bioconda star -y
# First create genome index
star --runMode genomeGenerate \
--genomeDir reference_genome_index \
--genomeFastaFiles reference_genome.fasta \
--sjdbGTFfile annotation.gtf \
--sjdbOverhang 99