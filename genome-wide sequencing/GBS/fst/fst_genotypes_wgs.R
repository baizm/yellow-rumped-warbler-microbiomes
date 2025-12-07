setwd("~/Documents/projects/YRWA_GWAS_microbiome/fst_genotypes/fst_genotypes2/")
library(ggplot2)

#gemmaWGS
#load genotype data for the significant SNPs
gtw<-read.table('gemmaWGS_genotypes_WGSdataset.txt')
idsw<-read.table('gemma_idsWGS_bgzip_ids.txt'); idsw<-substr(idsw$V1,1,7)
colnames(gtw)<-c('CHROM','POS','REF','ALT',idsw)

#reshape gt3 to make plotting easier
genosw<-gtw[5:ncol(gtw)]
genosw<-data.frame(t(genosw))
sites2<-paste(gtw$CHROM,gtw$POS,sep=':'); colnames(genosw)<-sites2
alpha_r2<-readRDS("~/Documents/projects/YRWA_GWAS_microbiome/qiime/phyloseq/alpha_r.rds") #read in alpha data
alpha_r2<-alpha_r2[which(rownames(alpha_r2) %in% rownames(genosw)),] 
identical(rownames(genosw), rownames(alpha_r2)) #TRUE
genosw$Chao1<-alpha_r2$Chao1; genosw$Shannon<-alpha_r2$Shannon; genosw$PD<-alpha_r2$PD

#turn the | into / so phased genotypes are grouped with unphased and fix NAs (missing genotypes)
genosw.1<-data.frame(lapply(genosw, function(x) {
  gsub("\\./.", NA, x)
}))

genosw.1$Chao1<-as.numeric(genosw.1$Chao1)
genosw.1$Shannon<-as.numeric(genosw.1$Shannon)
genosw.1$PD<-as.numeric(genosw.1$PD)

#RASGEF1C (alpha)
png('wgs_RASGEF1C.png', width=3.5, height=3, units='in', res=1200)
ggplot(subset(genosw.1, !is.na(chr13.7721206)), aes(x=chr13.7721206, y=Chao1)) + 
  geom_boxplot(outliers=F) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='gray60') +
  theme_bw() + theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 14))
dev.off()

# read in pheno data to plot pcoa 1+2
pheno<-read.table("~/Documents/projects/YRWA_GWAS_microbiome/gemma/gemmaWGS/run_ccr/in_file_gemmaWGS.fam")
genosw.1$pcoa1<-pheno$V19; genosw.1$pcoa2<-pheno$V20

#XDH
png('wgs_XDH.png', width=3.5, height=3, units='in', res=1200)
ggplot(subset(genosw.1, !is.na(chr3.6484377)), aes(x=chr3.6484377, y=pcoa2)) + 
  geom_boxplot(outliers=F) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='gray60') +
  theme_bw() + theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 14))
dev.off()

#the top beta1
png('wgs_chr11.15190314.png', width=3.5, height=3, units='in', res=1200)
ggplot(subset(genosw.1, !is.na(chr11.15190314)),aes(x=chr11.15190314, y=pcoa1)) + 
  geom_boxplot(outliers=F) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='gray60') +
  theme_bw() + theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 14)) #
dev.off()

#the top beta2
png('wgs_ELF1.png', width=3.5, height=3, units='in', res=1200)
ggplot(subset(genosw.1, !is.na(chr1.55000075)),aes(x=chr1.55000075, y=pcoa2)) + 
  geom_boxplot(outliers=F) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='gray60') +
  theme_bw() + theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 14)) #
dev.off()
