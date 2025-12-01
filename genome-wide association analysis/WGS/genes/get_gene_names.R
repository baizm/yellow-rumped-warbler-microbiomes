setwd("~/Documents/projects/YRWA_GWAS_microbiome/gemma/genes/genes_WGS/blastn/")
library(stringr)
        
load('get_gene_names.RData')

#read in TaeGut gtf file with gene names and ensemble ids to be matched
g<-read.table('~/Documents/projects/TaeGut_genome/Taeniopygia_guttata.bTaeGut1_v1.p.113.gtf.gz',skip=4,fill=T,sep=';')

g2<-data.frame(type=str_split_i(g$V1, "\t", 3),
               geneid=str_split_i(g$V1, " ", 2),
               trsid=g$V3,
               gene_name=g$V5)
g3<-g2[grep('ENSTGUT', g2$trsid),] #pull out only rows with transcript ids
g4<-g3[grep('gene_name', g3$gene_name),] #pull out only rows with gene_name (will leave others: "gene_source ensembl"--aka novel predicted genes w/o name, e.g. genes[which(genes$transcriptid=='ENSTGUT00000007346'),])
g5<-g4[,2:4]
g6<-g5[!duplicated(g5), ]

#-----finally, use this for searching!
genes<-data.frame(geneid=g6$geneid,
                  transcriptid=str_split_i(g6$trsid," ", 3),
                  gene=str_split_i(g6$gene_name," ", 3))

##read in top blast hits
#alpha
a<-read.table('taegut_genes_top_a.txt', header=T)
a$transcriptid<-substr(a$sseqid,1,18)

genes_a<-merge(x = a, y = genes, by = "transcriptid", all.x = TRUE) #includes transcripts w/o gene_name as NA
genes_a<-genes_a[!duplicated(genes_a$geneid), ] #get rid of genes with multiple transcripts
genes_a<-genes_a[!is.na(genes_a$gene), ] #get rid of transcripts w/NA gene_name

write.table(genes_a, 'candidate_genes_alpha.txt', sep='\t', row.names=F, col.names=F, quote=F)

#beta PCOA1
b1<-read.table('taegut_genes_top_b1.txt', header=T)
b1$transcriptid<-substr(b1$sseqid,1,18)

genes_b1<-merge(x = b1, y = genes, by = "transcriptid", all.x = TRUE) #includes transcripts w/o gene_name as NA
genes_b1<-genes_b1[!duplicated(genes_b1$geneid), ] #get rid of genes with multiple transcripts
genes_b1<-genes_b1[!is.na(genes_b1$gene), ] #get rid of transcripts w/NA gene_name

write.table(genes_b1, 'candidate_genes_beta1.txt', sep='\t', row.names=F, col.names=F, quote=F)

#beta PCOA2
b2<-read.table('taegut_genes_top_b2.txt', header=T)
b2$transcriptid<-substr(b2$sseqid,1,18)

genes_b2<-merge(x = b2, y = genes, by = "transcriptid", all.x = TRUE) #includes transcripts w/o gene_name as NA
genes_b2<-genes_b2[!duplicated(genes_b2$geneid), ] #get rid of genes with multiple transcripts
genes_b2<-genes_b2[!is.na(genes_b2$gene), ] #get rid of transcripts w/NA gene_name

write.table(genes_b2, 'candidate_genes_beta2.txt', sep='\t', row.names=F, col.names=F, quote=F)


#---make dataframes with gemma snps and warbler gff ids, too-----
alpha_results<-readRDS('~/Documents/projects/YRWA_GWAS_microbiome/gemma/genes/genes_WGS/alpha_results.RDS')
#add zebra finch ensemble gene ID and the gene name, matching the query sequence with warbler gff ID
alpha_results$geneid<-genes_a$geneid[match(alpha_results$warb_gff_geneID, genes_a$qseqid, )]
alpha_results$gene<-genes_a$gene[match(alpha_results$warb_gff_geneID, genes_a$qseqid, )]
saveRDS(alpha_results, 'alpha_results2.RDS') 

beta1_results<-readRDS('~/Documents/projects/YRWA_GWAS_microbiome/gemma/genes/genes_WGS/beta1_results.RDS')
#add zebra finch ensemble gene ID and the gene name, matching the query sequence with warbler gff ID
beta1_results$geneid<-genes_b1$geneid[match(beta1_results$warb_gff_geneID, genes_b1$qseqid, )]
beta1_results$gene<-genes_b1$gene[match(beta1_results$warb_gff_geneID, genes_b1$qseqid, )]
saveRDS(beta1_results, 'beta1_results2.RDS') 

beta2_results<-readRDS('~/Documents/projects/YRWA_GWAS_microbiome/gemma/genes/genes_WGS/beta2_results.RDS')
#add zebra finch ensemble gene ID and the gene name, matching the query sequence with warbler gff ID
beta2_results$geneid<-genes_b2$geneid[match(beta2_results$warb_gff_geneID, genes_b2$qseqid, )]
beta2_results$gene<-genes_b2$gene[match(beta2_results$warb_gff_geneID, genes_b2$qseqid, )]
saveRDS(beta2_results, 'beta2_results2.RDS') 

##---read in top blast hits for all genes---
all<-read.table('taegut_genes_top_all.txt', header=T)
all$transcriptid<-substr(all$sseqid,1,18) #cut off decimal

genes_all<-merge(x = all, y = genes, by = "transcriptid", all.x = TRUE) #includes transcripts w/o gene_name as NA
genes_all<-genes_all[!duplicated(genes_all$geneid), ] #get rid of genes with multiple transcripts
genes_all<-genes_all[!is.na(genes_all$gene), ] #get rid of transcripts w/NA gene_name

write.table(genes_all, 'candidate_genes_reference_list_gemmaWGS.txt', sep='\t', row.names=F, col.names=F, quote=F)
length(unique(genes_all$gene)) 



save.image('get_gene_names.RData')
