# run pear for merging paired end reads and alignment

/rhome/dpier004/pear/bin/pear -f indiv_R1_001.fastq.gz -r indiv_R2_001.fastq.gz -k -o indiv -j 24
rm indiv.discarded.fastq

#bwa paired reads
bwa mem -t 24 ~/shared/bwaDb/mywagenomev2.1.fa indiv.unassembled.forward.fastq indiv.unassembled.reverse.fastq | \
samtools fixmate -@ 2 -m - - | samtools sort -@ 2 -o - | samtools markdup -@ 2 -rs - indiv.pair.bam

rm indiv.unassembled.*.fastq

#bwa single reads
bwa mem -t 24 ~/shared/bwaDb/mywagenomev2.1.fa indiv.assembled.fastq | \
samtools sort -@ 2 - | samtools markdup -@ 2 -rs - indiv.single.bam

rm indiv.assembled.fastq

#merge the two bams
samtools merge indiv.all.bam indiv.pair.bam indiv.single.bam -@ 24

rm indiv.pair.bam indiv.single.bam

#index the final bam
samtools index indiv.all.bam
