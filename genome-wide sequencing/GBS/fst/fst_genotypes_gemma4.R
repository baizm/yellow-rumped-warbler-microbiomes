setwd("~/Documents/projects/YRWA_GWAS_microbiome/fst_genotypes/fst_gemma4/")
library(ggplot2)

load('fst_genotypes_gemma4.RData')

#--- FST and outliers ---
#load fst information
fst<-read.table('gemma_ids_no_thin_fst.txt.weir.fst', header=T) 
mean(fst$WEIR_AND_COCKERHAM_FST,na.rm = T)

#make vector of outlier fst snps 
q<-quantile(fst$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = T) 

outliers<-fst[which(fst$WEIR_AND_COCKERHAM_FST > q),]
outliers2<-paste(outliers$CHROM,outliers$POS,sep=':') #outliers in chr:pos format

#load genotype information gemma4
gt3<-read.table('gemma4_genotypes.txt')
ids3<-read.table('gemma_ids_no_thin_bgzip_ids.txt'); ids<-substr(ids3$V1,1,7)
colnames(gt3)<-c('CHROM','POS','REF','ALT',ids)
sig_snps<-paste(gt3$CHROM,gt3$POS,sep=':') #sig_snps in chr:pos format

#which gemma snps are FST outliers? 
outliers$gemma4_sig <- outliers2 %in% sig_snps
outliers[which(outliers$gemma4_sig==T),] #chr3:81754635


#----look at genotype ~ alpha for significant SNPs.
#reshape gt3 to make plotting easier
genos3<-gt3[5:ncol(gt3)]
genos3<-data.frame(t(genos3))
sites<-paste(gt3$CHROM,gt3$POS,sep=':'); colnames(genos3)<-sites

alpha_r<-readRDS("~/Documents/projects/YRWA_GWAS_microbiome/qiime/phyloseq/alpha_r.rds") #read in alpha data
alpha_r<-alpha_r[which(rownames(alpha_r) %in% rownames(genos3)),] 
identical(rownames(genos3), rownames(alpha_r)) #TRUE
genos3$Chao1<-alpha_r$Chao1; genos3$Shannon<-alpha_r$Shannon; genos3$PD<-alpha_r$PD

#turn the | into / so phased genotypes are grouped with unphased and fix NAs (missing genotypes)
genos3.1<-data.frame(lapply(genos3, function(x) {
  gsub("\\|", "/", x)
}))
genos3.1<-data.frame(lapply(genos3.1, function(x) {
  gsub("\\./.", NA, x)
}))

genos3.1$Chao1<-as.numeric(genos3.1$Chao1)
genos3.1$Shannon<-as.numeric(genos3.1$Shannon)
genos3.1$PD<-as.numeric(genos3.1$PD)

# to check which allele is which C=myrtle (compare to admixture classification for those individuals)
table(genos3.1$chr3.81754635)
data.frame(row.names(genos3),genos3$`chr3:81754635`) 
#read in myr and aud ids
myr_ids<-read.table('myr_ids2.txt'); myr_ids$V1<-substr(myr_ids$V1,1,7)
aud_ids<-read.table('aud_ids2.txt'); aud_ids$V1<-substr(aud_ids$V1,1,7)
chr3.81754635 <- data.frame(id=row.names(genos3),gt.chr3.81754635=genos3$`chr3:81754635`)
chr3.81754635$species <- ifelse(chr3.81754635$id %in% myr_ids$V1,'myr',
                                ifelse(chr3.81754635$id %in% aud_ids$V1, 'aud', 'hyb'))
#find aud allele:
chr3.81754635[which(chr3.81754635$species=='aud'),]
chr3.81754635[which(chr3.81754635$species=='myr'),]

table(chr3.81754635$gt.chr3.81754635[which(chr3.81754635$species=='aud')])
table(chr3.81754635$gt.chr3.81754635[which(chr3.81754635$species=='myr')])
table(chr3.81754635$gt.chr3.81754635[which(chr3.81754635$species=='hyb')])


#write table of fst for significant snps in gemma4
fst_sig_snps<-fst; fst_sig_snps$chr.pos<-paste(fst$CHROM,fst$POS,sep=":")
fst_sig_snps<-fst_sig_snps[which(fst_sig_snps$chr.pos %in% sig_snps),]
write.csv(fst_sig_snps, 'fst_sig_snps.csv', quote=F,col.names = T, row.names = F, na = "NA")


#read in phenotype data
#note: chao1 values plotted are not log-transformed...they are the original
pheno<-read.table("~/Documents/projects/YRWA_GWAS_microbiome/gemma/gemma4/in_file_gemma.fam")
genos3.1$pcoa1<-pheno$V19; genos3.1$pcoa2<-pheno$V20 
pheno2<-read.table("~/Documents/projects/YRWA_GWAS_microbiome/gemma/gemma3/gemmaASVs/asvs_for_gemma3_pheno.txt")
genos3.1$beij<-pheno2$V1; genos3.1$micro<-pheno2$V2

#plot abundance of ASVs by genotype

#beta1
png('figs/gemma4_chrz.61657428.png', width=3.5, height=3, units='in', res=1200)
ggplot(subset(genos3.1, !is.na(chrz.61657428)), aes(x=chrz.61657428, y=pcoa1)) + 
  geom_boxplot(outliers=F) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='gray60') +
  theme_bw() + theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 14)) #NIPBL
dev.off()


#the one FST outlier
png('figs/gemma4_chr3.81754635.png', width=3.5, height=3, units='in', res=1200)
ggplot(subset(genos3.1, !is.na(chr3.81754635)), aes(x=chr3.81754635, y=pcoa2)) + 
  geom_boxplot(outliers=F) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='gray60') +
  theme_bw() + theme(axis.text = element_text(size = 14)) + theme(axis.title = element_text(size = 14))#no gene
dev.off()


save.image('fst_genotypes_gemma4.RData')