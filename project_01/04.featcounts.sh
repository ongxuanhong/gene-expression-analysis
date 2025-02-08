# Create directory for counts
mkdir -p counts

# Run featureCounts
featureCounts -p -t gene -g gene_id \
featureCounts -p -t gene -g gene_id -s 1 \
  -a gencode.v47lift37.annotation.gtf \
  -o counts/counts.txt \
  alignments/Cont1.sorted.bam \
  alignments/Zn2.sorted.bam  