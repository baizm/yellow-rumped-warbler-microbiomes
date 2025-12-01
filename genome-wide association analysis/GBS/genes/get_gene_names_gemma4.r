setwd("~/Documents/projects/YRWA_GWAS_microbiome/gemma/genes/genes_gemma4/blastn/")
library(stringr)

load('get_gene_names_gemma4.RData')

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
#beta PCOA1
b1<-read.table('taegut_genes_top_b1.txt', header=T)
b1$transcriptid<-substr(b1$sseqid,1,18)

genes_b1<-merge(x = b1, y = genes, by = "transcriptid", all.x = TRUE) #includes transcripts w/o gene_name as NA
genes_b1<-genes_b1[!duplicated(genes_b1$geneid), ] #get rid of genes with multiple transcripts
genes_b1<-genes_b1[!is.na(genes_b1$gene), ] #get rid of transcripts w/NA gene_name

write.table(genes_b1, 'candidate_genes_beta1_gemma4.txt', sep='\t', row.names=F, col.names=F, quote=F)

#beij
beij<-read.table('taegut_genes_top_beij.txt', header=T)
beij$transcriptid<-substr(beij$sseqid,1,18)

genes_beij<-merge(x = beij, y = genes, by = "transcriptid", all.x = TRUE) #includes transcripts w/o gene_name as NA
genes_beij<-genes_beij[!duplicated(genes_beij$geneid), ] #get rid of genes with multiple transcripts
genes_beij<-genes_beij[!is.na(genes_beij$gene), ] #get rid of transcripts w/NA gene_name

write.table(genes_beij, 'candidate_genes_beij_gemma4.txt', sep='\t', row.names=F, col.names=F, quote=F)

#micro
micro<-read.table('taegut_genes_top_micro.txt', header=T)
micro$transcriptid<-substr(micro$sseqid,1,18)

genes_micro<-merge(x = micro, y = genes, by = "transcriptid", all.x = TRUE) #includes transcripts w/o gene_name as NA
genes_micro<-genes_micro[!duplicated(genes_micro$geneid), ] #get rid of genes with multiple transcripts
genes_micro<-genes_micro[!is.na(genes_micro$gene), ] #get rid of transcripts w/NA gene_name

write.table(genes_micro, 'candidate_genes_micro_gemma4.txt', sep='\t', row.names=F, col.names=F, quote=F)

#all genes for go enrichment
all<-read.table('taegut_genes_top_all.txt', header=T)
all$transcriptid<-substr(all$sseqid,1,18)

genes_all<-merge(x = all, y = genes, by = "transcriptid", all.x = TRUE) #includes transcripts w/o gene_name as NA
genes_all<-genes_all[!duplicated(genes_all$geneid), ] #get rid of genes with multiple transcripts
genes_all<-genes_all[!is.na(genes_all$gene), ] #get rid of transcripts w/NA gene_name

write.table(genes_all, 'candidate_genes_all_gemma4.txt', sep='\t', row.names=F, col.names=F, quote=F)

#now write file for panther gene reference list
full4<-genes_all$gene; length(full4) #8895 genes
full4<-unique(full4); length(full4) #8895 genes; none duplicated, test this one for the reference set

write.table(full4, './GO_enrichment/candidate_genes_all_gemma4_reference_genes.txt',row.names = F, quote=F, col.names = F)

save.image('get_gene_names_gemma4.RData')
