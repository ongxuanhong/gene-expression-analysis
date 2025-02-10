# Create directory for counts
mkdir -p counts

# Run featureCounts
featureCounts -p -t exon -g gene_id \
  -a hisat2_index/gencode.v19.annotation.gtf \
  -o counts/counts.txt \
  alignments/Cont1.sorted.bam \
  alignments/Zn2.sorted.bam  