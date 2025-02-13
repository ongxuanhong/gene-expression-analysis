# Create directory for counts
mkdir -p counts

featureCounts -p -t exon -g gene_id \
    -a hisat2_index/gencode.v19.annotation.gtf \
    -o counts/gene_counts.tsv \
    bam_files/*.sorted.bam  