setwd("~/Documents/projects/YRWA_GWAS_microbiome/gemma/genes/genes_gemma4/")
load('get_genes_gemma4.RData')

#load sig_snps from gemma4
sig_snps<-read.csv('../../gemma4/sig_snps.csv') #written from man_plots_gemma4.R

#---significant gemma snps----
regions_a<-sig_snps[which(sig_snps$pheno=='alpha'),c('chr','ps')]
colnames(regions_a)<-c('chr','pos')
table(regions_a$chr)

#beta PCOA1 regions
regions_b1<-sig_snps[which(sig_snps$pheno=='beta1'),c('chr','ps')]
colnames(regions_b1)<-c('chr','pos')

#beta PCOA2 regions
regions_b2<-sig_snps[which(sig_snps$pheno=='beta2'),c('chr','ps')]
colnames(regions_b2)<-c('chr','pos')

#beij regions
regions_beij<-sig_snps[which(sig_snps$pheno=='beij'),c('chr','ps')]
colnames(regions_beij)<-c('chr','pos')

#micro regions
regions_micro<-sig_snps[which(sig_snps$pheno=='micro'),c('chr','ps')]
colnames(regions_micro)<-c('chr','pos')

#----define 20KB buffer window around significant snps/genes actually 
buf<-20000

#---read in the gene annotation file
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

#beij
regions_beij$match= sapply( 1:nrow(regions_beij) , 
                          function(x)   
                            any(  regions_beij[x, 'chr']==gff[, 'chr'] &
                                    regions_beij[x , 'pos'] <= gff[ , 'end2'] & 
                                    regions_beij[x , 'pos'] >= gff[ , 'start2'] ))
table(regions_beij$match) 

#micro
regions_micro$match= sapply( 1:nrow(regions_micro) , 
                          function(x)   
                            any(  regions_micro[x, 'chr']==gff[, 'chr'] &
                                    regions_micro[x , 'pos'] <= gff[ , 'end2'] & 
                                    regions_micro[x , 'pos'] >= gff[ , 'start2'] ))
table(regions_micro$match) 

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
gff[which(gff$match_b1==T),]

#beta PCOA2
gff$match_b2= sapply( 1:nrow(gff) , 
                      function(x)   
                        any(  gff[x, 'chr']==regions_b2[, 'chr'] &
                                gff[x , 'end2'] >= regions_b2[ , 'pos'] & 
                                gff[x , 'start2'] <= regions_b2[ , 'pos'] ))
table(gff$match_b2) 

#beij
gff$match_beij= sapply( 1:nrow(gff) , 
                      function(x)   
                        any(  gff[x, 'chr']==regions_beij[, 'chr'] &
                                gff[x , 'end2'] >= regions_beij[ , 'pos'] & 
                                gff[x , 'start2'] <= regions_beij[ , 'pos'] ))
table(gff$match_beij) 

#micro
gff$match_micro= sapply( 1:nrow(gff) , 
                      function(x)   
                        any(  gff[x, 'chr']==regions_micro[, 'chr'] &
                                gff[x , 'end2'] >= regions_micro[ , 'pos'] & 
                                gff[x , 'start2'] <= regions_micro[ , 'pos'] ))
table(gff$match_micro) 


dim(gff[which(gff$match_a==T & gff$match_b1==T),])
gff[which(gff$match_a==T & gff$match_b1==T),] #0 genes sig for a and b1
gff[which(gff$match_a==T & gff$match_b2==T),] #0 genes sig for a and b2
gff[which(gff$match_a==T & gff$match_beij==T),] #0 genes sig for a and beij
gff[which(gff$match_a==T & gff$match_micro==T),] #0 genes sig for a and micro

gff[which(gff$match_b1==T & gff$match_b2==T),] #0 genes sig for a and b1
gff[which(gff$match_b1==T & gff$match_beij==T),] #0 genes sig for a and beij
gff[which(gff$match_b1==T & gff$match_micro==T),] #0 genes sig for a and micro

gff[which(gff$match_b2==T & gff$match_beij==T),] #0 genes sig for b2 and beij
gff[which(gff$match_b2==T & gff$match_micro==T),] #0 genes sig for b2 and micro

gff[which(gff$match_beij==T & gff$match_micro==T),] #0 genes sig for beij and micro

##-------pull out IDS for unique genes----
#alpha genes
gffa<-gff[which(gff$match_a==T),]
ids_a<-data.frame(do.call('rbind', strsplit(as.character(gffa$attr),';',fixed=TRUE))) #split attr column
ids_a2<-substring(ids_a$X1,4)
#no genes

#b1 genes
gffb1<-gff[which(gff$match_b1==T),]
ids_b1<-data.frame(do.call('rbind', strsplit(as.character(gffb1$attr),';',fixed=TRUE))) #split attr column
ids_b1_2<-substring(ids_b1$X1,4)
write.table(ids_b1_2, 'gene_ids_20kb_b1.txt', quote=F, col.names = F,row.names=F)

#b2 genes
gffb2<-gff[which(gff$match_b2==T),]
ids_b2<-data.frame(do.call('rbind', strsplit(as.character(gffb2$attr),';',fixed=TRUE))) #split attr column
ids_b2_2<-substring(ids_b2$X1,4)
#no genes

#beij genes
gffbeij<-gff[which(gff$match_beij==T),]
ids_beij<-data.frame(do.call('rbind', strsplit(as.character(gffbeij$attr),';',fixed=TRUE))) #split attr column
ids_beij_2<-substring(ids_beij$X1,4)
write.table(ids_beij_2, 'gene_ids_20kb_beij.txt', quote=F, col.names = F,row.names=F)

#micro genes
gffmicro<-gff[which(gff$match_micro==T),]
ids_micro<-data.frame(do.call('rbind', strsplit(as.character(gffmicro$attr),';',fixed=TRUE))) #split attr column
ids_micro_2<-substring(ids_micro$X1,4)
write.table(ids_micro_2, 'gene_ids_20kb_micro.txt', quote=F, col.names = F,row.names=F)

##next: 
#used seqkit module on CCR:
#seqkit grep -f gene_ids_20kb_a.txt /projects/academic/mbaiz/mywa_genome/annotation/mywagenomev2.1.all.maker.transcripts.fasta > gene_seqs_20kb_a.fasta
#pull out sequences, blasted against zebra finch cds (used blastdb from cds downloaded from ensembl for WGS dataset)
#used get_top_hit.py to parse the blastn results

#---Get All gemma snps for reference GO enrichment (to comapre with significant SNPs)----
#all snps tested 
regions_all<-read.table('../../gemma4/output/gemma_out4_a.assoc.txt.gz', header=T) #same number of snps for b1 and b2 as well, can just use this file for all
table(regions_all$chr)
#note: don't have to get rid of scaffolds to match gff and gemma3 snps
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

#used seqkit module on CCR:
#seqkit grep -f gene_ids_20kb_all.txt /projects/academic/mbaiz/mywa_genome/annotation/mywagenomev2.1.all.maker.transcripts.fasta > gene_seqs_20kb_all.fasta
#pull out sequences, blasted against zebra finch cds (used blastdb from cds downloaded from ensembl for WGS dataset)
#used get_top_hit.py to parse the blastn results


#get rest of info for supp table
#alpha - none
#beta2 - none

#beta1
table(regions_b1$match) #4 TRUE SNPs, 4 FALSE SNPs
gff[which(gff$match_b1==T),] #6 predicted genes

#beij
table(regions_beij$match) #2 TRUE SNPs, 10 FALSE SNPs
gff[which(gff$match_beij==T),] #2 predicted genes

#micro
table(regions_micro$match) #1 TRUE SNPs, 2 FALSE SNPs
gff[which(gff$match_micro==T),] #1 predicted genes

gff_for_supp<-rbind(gff[which(gff$match_b1==T),], gff[which(gff$match_beij==T),], gff[which(gff$match_micro==T),])
write.csv(gff_for_supp, 'gff_for_supp.csv', row.names = F, quote=F) #manually pasted into Supplementary_table file

save.image('get_genes_gemma4.RData')

