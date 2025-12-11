# run bcftools mpileup on one portion of the genome. Requires a file bamlist.txt with paths to all bam files.

bcftools mpileup -Ou -q 20 -a DP -d 500 -R beds/chunk.bed \
-f /bigdata/brelsfordlab/shared/bwaDb/mywagenomev2.1.fa -b bamlist.txt | \
bcftools call -Oz -mv -f GQ -o output/chunk.vcf.gz
