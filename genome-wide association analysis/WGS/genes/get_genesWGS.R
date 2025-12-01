setwd("~/Documents/projects/YRWA_GWAS_microbiome/gemma/genes/genes_WGS/")
load('get_genesWGS.RData')

#alpha regions 
regions_a<-read.table('../../gemmaWGS/output/snps_a_gemmaWGS.txt')
colnames(regions_a)<-c('chr','pos')
table(regions_a$chr)

#beta PCOA1 regions 
regions_b1<-read.table('../../gemmaWGS/output/snps_b1_gemmaWGS.txt')
colnames(regions_b1)<-c('chr','pos')

#beta PCOA2 regions 
regions_b2<-read.table('../../gemmaWGS/output/snps_b2_gemmaWGS.txt')
colnames(regions_b2)<-c('chr','pos')

#define 20KB buffer window around significant snps/genes actually 
buf<-20000

#read in the gene annotation file
gff<-read.table('../mywagenomev2.1.all.noseq.gff.gz', sep='\t',)
colnames(gff)<-c('chr','source','type','start','end','score','strand','phase','attr')
gff<-gff[which(gff$type=='mRNA'),] #keep protein coding only
gff<-gff[grep('chr', gff$chr),] #get rid of scaffolds and mito
unique(gff$chr[order(gff$chr)])
#add 20KB buffer
gff$start2<-gff$start-buf
gff$end2<-gff$end+buf

##---find how many significant SNPs are near genes-----------
#alpha
regions_a$match= sapply( 1:nrow(regions_a) , 
                   function(x)   
                     any(  regions_a[x, 'chr']==gff[, 'chr'] &
                             regions_a[x , 'pos'] <= gff[ , 'end2'] & 
                             regions_a[x , 'pos'] >= gff[ , 'start2'] ))
table(regions_a$match) 

#beta PCOA1
regions_b1$match= sapply( 1:nrow(regions_b1) , 
                         function(x)   
                           any(  regions_b1[x, 'chr']==gff[, 'chr'] &
                                   regions_b1[x , 'pos'] <= gff[ , 'end2'] & 
                                   regions_b1[x , 'pos'] >= gff[ , 'start2'] ))
table(regions_b1$match) 

#beta PCOA2
regions_b2$match= sapply( 1:nrow(regions_b2) , 
                          function(x)   
                            any(  regions_b2[x, 'chr']==gff[, 'chr'] &
                                    regions_b2[x , 'pos'] <= gff[ , 'end2'] & 
                                    regions_b2[x , 'pos'] >= gff[ , 'start2'] ))
table(regions_b2$match) 


###---find gff GENES whose 20KB buffer overlaps a SNP-------
#alpha
gff$match_a= sapply( 1:nrow(gff) , 
                         function(x)   
                           any(  gff[x, 'chr']==regions_a[, 'chr'] &
                                   gff[x , 'end2'] >= regions_a[ , 'pos'] & 
                                   gff[x , 'start2'] <= regions_a[ , 'pos'] ))
table(gff$match_a) 

#beta PCOA1
gff$match_b1= sapply( 1:nrow(gff) , 
                     function(x)   
                       any(  gff[x, 'chr']==regions_b1[, 'chr'] &
                               gff[x , 'end2'] >= regions_b1[ , 'pos'] & 
                               gff[x , 'start2'] <= regions_b1[ , 'pos'] ))
table(gff$match_b1) 

#beta PCOA2
gff$match_b2= sapply( 1:nrow(gff) , 
                      function(x)   
                        any(  gff[x, 'chr']==regions_b2[, 'chr'] &
                                gff[x , 'end2'] >= regions_b2[ , 'pos'] & 
                                gff[x , 'start2'] <= regions_b2[ , 'pos'] ))
table(gff$match_b2) 

dim(gff[which(gff$match_a==T & gff$match_b1==T),])
gff[which(gff$match_a==T & gff$match_b1==T),]

##-------pull out IDS for unique genes----
#alpha genes
gffa<-gff[which(gff$match_a==T),]
ids_a<-data.frame(do.call('rbind', strsplit(as.character(gffa$attr),';',fixed=TRUE))) #split attr column
ids_a2<-substring(ids_a$X1,4)
write.table(ids_a2, 'gene_ids_20kb_a.txt', quote=F, col.names = F,row.names=F)

#b1 genes
gffb1<-gff[which(gff$match_b1==T),]
ids_b1<-data.frame(do.call('rbind', strsplit(as.character(gffb1$attr),';',fixed=TRUE))) #split attr column
ids_b1_2<-substring(ids_b1$X1,4)
write.table(ids_b1_2, 'gene_ids_20kb_b1.txt', quote=F, col.names = F,row.names=F)

#b2 genes
gffb2<-gff[which(gff$match_b2==T),]
ids_b2<-data.frame(do.call('rbind', strsplit(as.character(gffb2$attr),';',fixed=TRUE))) #split attr column
ids_b2_2<-substring(ids_b2$X1,4)
write.table(ids_b2_2, 'gene_ids_20kb_b2.txt', quote=F, col.names = F,row.names=F)


##next:
#used seqkit module on CCR
#seqkit grep -f gene_ids_20kb_a.txt /projects/academic/mbaiz/mywa_genome/annotation/mywagenomev2.1.all.maker.transcripts.fasta > gene_seqs_20kb_a.fasta
#pull out sequences, blated against zebra finch cds (made blastdb from cds downloaded from ensembl)

#---associate gff genes with original gemma3 significant SNPs---
#alpha
alpha_results<-gffa #copy gffa data.frame
alpha_results$warb_gff_geneID <- ids_a2 #add predicted gene ID from annotated genome
#add column of matching gemma sig snps (returns list of row indexes for matching regions_a)
alpha_results$gemma_snp= sapply( 1:nrow(gffa) , 
                                 function(x)   
                                   which(  gffa[x, 'chr']==regions_a[, 'chr'] &
                                             gffa[x , 'end2'] >= regions_a[ , 'pos'] & 
                                             gffa[x , 'start2'] <= regions_a[ , 'pos'] ))
alpha_results<-tidyr::unnest(alpha_results, gemma_snp) #need to split rows where multiple snps match the same gff gene
alpha_results$gemma_chr<-regions_a$chr[alpha_results$gemma_snp] #add gemma chr by row index
alpha_results$gemma_pos<-regions_a$pos[alpha_results$gemma_snp] #add gemma pos by row index
saveRDS(alpha_results, 'alpha_results.RDS')

#beta pcoa1
beta1_results<-gffb1 #copy gffb1 data.frame
beta1_results$warb_gff_geneID <- ids_b1_2 #add predicted gene ID from annotated genome
#add column of matching gemma sig snps 
beta1_results$gemma_snp= sapply( 1:nrow(gffb1) , 
                                 function(x)   
                                   which(  gffb1[x, 'chr']==regions_b1[, 'chr'] &
                                             gffb1[x , 'end2'] >= regions_b1[ , 'pos'] & 
                                             gffb1[x , 'start2'] <= regions_b1[ , 'pos'] ))
beta1_results<-tidyr::unnest(beta1_results, gemma_snp) #need to split rows where multiple snps match the same gff gene
beta1_results$gemma_chr<-regions_b1$chr[beta1_results$gemma_snp] #add gemma chr by row index
beta1_results$gemma_pos<-regions_b1$pos[beta1_results$gemma_snp] #add gemma pos by row index
saveRDS(beta1_results, 'beta1_results.RDS')

#beta pcoa2
beta2_results<-gffb2 #copy gffb2 data.frame
beta2_results$warb_gff_geneID <- ids_b2_2 #add predicted gene ID from annotated genome
#add column of matching gemma sig snps 
beta2_results$gemma_snp= sapply( 1:nrow(gffb2) , 
                                 function(x)   
                                   which(  gffb2[x, 'chr']==regions_b2[, 'chr'] &
                                             gffb2[x , 'end2'] >= regions_b2[ , 'pos'] & 
                                             gffb2[x , 'start2'] <= regions_b2[ , 'pos'] ))
beta2_results<-tidyr::unnest(beta2_results, gemma_snp) #need to split rows where multiple snps match the same gff gene
beta2_results$gemma_chr<-regions_b2$chr[beta2_results$gemma_snp] #add gemma chr by row index
beta2_results$gemma_pos<-regions_b2$pos[beta2_results$gemma_snp] #add gemma pos by row index
saveRDS(beta2_results, 'beta2_results.RDS')


#---Get All gemma snps for reference GO enrichment (to compare with significant SNPs)----
#all snps tested 
regions_all<-read.table('../../gemmaWGS/output/gemma_outWGS_a.assoc.txt.gz', header=T) #same number of snps for b1 and b2 as well, can just use this file for all
table(regions_all$chr)
#note: don't have to get rid of scaffolds to match gff and gemma snps
#fix chr names
regions_all$chr<-gsub("^chr", "", regions_all$chr)
regions_all$chr<-paste('chr',regions_all$chr, sep='')

##---find how many total SNPs are near genes-----------
regions_all$match= sapply( 1:nrow(regions_all) , 
                           function(x)   
                             any(  regions_all[x, 'chr']==gff[, 'chr'] &
                                     regions_all[x , 'ps'] <= gff[ , 'end2'] & 
                                     regions_all[x , 'ps'] >= gff[ , 'start2'] ))
table(regions_all$match) 

###---find gff GENES whose 20KB buffer overlaps a SNP-------
gff$match_all= sapply( 1:nrow(gff) , 
                       function(x)   
                         any(  gff[x, 'chr']==regions_all[, 'chr'] &
                                 gff[x , 'end2'] >= regions_all[ , 'ps'] & 
                                 gff[x , 'start2'] <= regions_all[ , 'ps'] ))
table(gff$match_all) 

##-------pull out IDS for unique genes----
#all genes
gffall<-gff[which(gff$match_all==T),]
ids_all<-data.frame(do.call('rbind', strsplit(as.character(gffall$attr),';',fixed=TRUE))) #split attr column
ids_all2<-substring(ids_all$X1,4)
write.table(ids_all2, 'gene_ids_20kb_all.txt', quote=F, col.names = F,row.names=F)

##next:
#used seqkit module on CCR:
#seqkit grep -f gene_ids_20kb_all.txt /projects/academic/mbaiz/mywa_genome/annotation/mywagenomev2.1.all.maker.transcripts.fasta > gene_seqs_20kb_all.fasta
#pull out sequences, blasted against zebra finch cds (used blastdb from cds downloaded from ensembl for WGS dataset)
#used get_top_hit.py to parse the blastn results


save.image('get_genesWGS.RData')
