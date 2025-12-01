# get fasta file of genes within 20KB of signfiicant SNPs

module load seqkit/2.7.0

seqkit grep -f gene_ids_20kb_b1.txt \
/projects/academic/mbaiz/mywa_genome/annotation/mywagenomev2.1.all.maker.transcripts.fasta > gene_seqs_20kb_b1.fasta

seqkit grep -f gene_ids_20kb_beij.txt \
/projects/academic/mbaiz/mywa_genome/annotation/mywagenomev2.1.all.maker.transcripts.fasta > gene_seqs_20kb_beij.fasta

seqkit grep -f gene_ids_20kb_micro.txt \
/projects/academic/mbaiz/mywa_genome/annotation/mywagenomev2.1.all.maker.transcripts.fasta > gene_seqs_20kb_micro.fasta
