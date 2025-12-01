##--script for plotting results of GEMMA analysis of GBS dataset----
setwd("~/Documents/projects/YRWA_GWAS_microbiome/gemma/gemma4/")

library(qqman)
library(qvalue)
library(dplyr)

load('man_plots_gemma4.RData')

###-----alpha----
a<-read.table('output/gemma_out4_a.assoc.txt.gz', header=T)

qa <- qvalue(p = a$p_lrt)
a$qvalue <- qa$qvalues
unique(a$chr) #1a,4a and z are all at the end
tail(a[which(a$chr=='chrz'),]) #1210435 is last chrz
a<-head(a, 1210435)
#re-order by chromosomes how I want
chrs_ord<-unique(a$chr); chrs_ord<-chrs_ord[c(1,29,2:4,30,5:28,31)]
a <- a %>% arrange(factor(chr, levels = chrs_ord))
#copy and manipulate further
a2<-mutate(a, chr2 = as.numeric(factor(chr, levels = chrs_ord)))
#make vector of sig snps to highlight
a2$rs2<-paste(rep('rs',length(a2$rs)),1:length(a2$rs),sep='') 
h_alpha<-a2$rs2[which(a2$p_lrt<1.140475e-06)] #threshold from permutations
#plot!
chrs_lab<-chrs_ord; chrs_lab[2]<-'1a'; chrs_lab[6]<-'4a'; chrs_lab[31]<-'Z' #make labels for plots
png('figs/man_alpha_gemma4.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(a2, chr = "chr2", bp = "ps", p = "p_lrt", snp = "rs2", chrlabs = chrs_lab, ylim=c(0,7),
          highlight=h_alpha, main='Chao 1',suggestiveline=F,
          cex=0.5, cex.axis=0.6)
abline(h=-log10(1.140475e-06), lty=2, col='palegreen3', lwd=0.8)
dev.off()

#significant SNPs
a2[which(a2$rs2 %in% h_alpha),]


##----beta PCOA1-----
b1<-read.table('output/gemma_out4_b1.assoc.txt.gz', header=T)

qb1 <- qvalue(p = b1$p_lrt)
b1$qvalue <- qb1$qvalues
unique(b1$chr) #1a,4a and z are all at the end
tail(b1[which(b1$chr=='chrz'),]) #223837 is last chrz

b1<-head(b1, 1210435)
#re-order by chromosomes how I want
b1 <- b1 %>% arrange(factor(chr, levels = chrs_ord))
#copy and manipulate further
b12<-mutate(b1, chr2 = as.numeric(factor(chr, levels = chrs_ord)))
#make vector of sig snps to highlight
b12$rs2<-paste(rep('rs',length(b12$rs)),1:length(b12$rs),sep='') 
h_beta1<-b12$rs2[which(b12$p_lrt<2.195969e-06)] #threshold from permutations
#plot!
png('figs/man_beta1_gemma4.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(b12, chr = "chr2", bp = "ps", p = "p_lrt", snp = "rs2", chrlabs = chrs_lab, ylim=c(0,7),
          highlight=h_beta1, main='UniFrac PCoA 1',suggestiveline=F,
          cex=0.5, cex.axis=0.6)
abline(h=-log10(2.195969e-06), lty=2, col='palegreen3', lwd=0.8)
dev.off()

#significant SNPs
b12[which(b12$rs2 %in% h_beta1),] 


##----beta PCOA2-----
b2<-read.table('output/gemma_out4_b2.assoc.txt.gz', header=T)

qb2 <- qvalue(p = b2$p_lrt)
b2$qvalue <- qb2$qvalues
unique(b2$chr) #1a,4a and z are all at the end
tail(b2[which(b2$chr=='chrz'),]) #1210435 is last chrz

b2<-head(b2, 1210435)
#re-order by chromosomes how I want
b2 <- b2 %>% arrange(factor(chr, levels = chrs_ord))
#copy and manipulate further
b22<-mutate(b2, chr2 = as.numeric(factor(chr, levels = chrs_ord)))
#make vector of sig snps to highlight
b22$rs2<-paste(rep('rs',length(b22$rs)),1:length(b22$rs),sep='') 
h_beta2<-b22$rs2[which(b22$p_lrt<2.829382e-06)] #threshold from permutations
#plot!
png('figs/man_beta2_gemma4.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(b22, chr = "chr2", bp = "ps", p = "p_lrt", snp = "rs2", chrlabs = chrs_lab, ylim=c(0,7.2),
          highlight=h_beta2, main='UniFrac PCoA 2',suggestiveline=F,
          cex=0.5, cex.axis=0.6)
abline(h=-log10(2.829382e-06), lty=2, col='palegreen3', lwd=0.8)
dev.off()

#significant SNPs
b22[which(b22$rs2 %in% h_beta2),] 


###-----Beijerinckiaceae ASV----
beij<-read.table('output/gemma_out4_beij.assoc.txt.gz', header=T)

qbeij <- qvalue(p = beij$p_lrt)
beij$qvalue <- qbeij$qvalues
unique(beij$chr) #1a,4a and z are all at the end
tail(beij[which(beij$chr=='chrz'),]) #1210435 is last chrz

beij<-head(beij, 1210435)
#re-order by chromosomes how I want
beij <- beij %>% arrange(factor(chr, levels = chrs_ord))
#copy and manipulate further
beij2<-mutate(beij, chr2 = as.numeric(factor(chr, levels = chrs_ord)))
#make vector of sig snps to highlight
beij2$rs2<-paste(rep('rs',length(beij2$rs)),1:length(beij2$rs),sep='') 
h_beij<-beij2$rs2[which(beij2$p_lrt<2.893159e-09)] #threshold from permutations

#plot!
png('figs/man_beij_gemma4.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(beij2, chr = "chr2", bp = "ps", p = "p_lrt", snp = "rs2", chrlabs = chrs_lab, ylim=c(0,13.5),
          highlight=h_beij, main='Beijerinckiaceae ASV',suggestiveline=F, genomewideline = F, cex=0.5, cex.axis=0.6)
abline(h=-log10(2.893159e-09), lty=2, col='palegreen3', lwd=0.8)
dev.off()

#significant SNPs
beij2[which(beij2$rs2 %in% h_beij),]


###-----micro----
micro<-read.table('output/gemma_out4_micro.assoc.txt.gz', header=T)

qmicro <- qvalue(p = micro$p_lrt)
micro$qvalue <- qmicro$qvalues
unique(micro$chr) #1a,4a and z are all at the end
tail(micro[which(micro$chr=='chrz'),]) #1210435 is last chrz

micro<-head(micro, 1210435)
#re-order by chromosomes how I want
micro <- micro %>% arrange(factor(chr, levels = chrs_ord))
#copy and manipulate further
micro2<-mutate(micro, chr2 = as.numeric(factor(chr, levels = chrs_ord)))
#make vector of sig snps to highlight
micro2$rs2<-paste(rep('rs',length(micro2$rs)),1:length(micro2$rs),sep='') 
h_micro<-micro2$rs2[which(micro2$p_lrt<8.776249e-07)] #threshold from permutations

#plot!
png('figs/man_micro_gemma4.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(micro2, chr = "chr2", bp = "ps", p = "p_lrt", snp = "rs2", chrlabs = chrs_lab, ylim=c(0,7),
          highlight=h_micro, main='Microbacteriaceae ASV',suggestiveline=F, genomewideline = F,
          cex=0.5, cex.axis=0.6)
abline(h=-log10(8.776249e-07), lty=2, col='palegreen3', lwd=0.8)
dev.off()

#significant SNPs
micro2[which(micro2$rs2 %in% h_micro),] 


#---build table of all significant SNPs---
sig_snps<-data.frame(rbind(a2[which(a2$rs2 %in% h_alpha),],
                           b12[which(b12$rs2 %in% h_beta1),],
                           b22[which(b22$rs2 %in% h_beta2),],
                           beij2[which(beij2$rs2 %in% h_beij),],
                           micro2[which(micro2$rs2 %in% h_micro),]))
sig_snps$pheno<-c(rep('alpha',3),rep('beta1',8),rep('beta2',12),rep('beij',12),rep('micro',3))

#order by position
sig_snps[with(sig_snps, order(chr, ps)), ]

#write sig_snps for get_genes_gemma4.r
write.csv(sig_snps, 'sig_snps.csv')

#write sig_snps for fst_gemma4 genotypes, add chr prefix
regions_gemma4<-sig_snps[,c('chr','ps')]
chrs_w_label<-ifelse(
  !startsWith(regions_gemma4$chr, 'chr'), # Check if the string already starts with the prefix
  paste0('chr', regions_gemma4$chr),     # If not, add the prefix
  regions_gemma4$chr                      # If it does, keep the original string
)

#check
regions_gemma4$chr; chrs_w_label #looks good, now paste new column
regions_gemma4$chr <- chrs_w_label
regions_gemma4

write.table(regions_gemma4, '../../fst_genotypes/fst_gemma4/regions_gemma4.txt', row.names = F, col.names = F, sep='\t', quote=F)

#--- make QQ plot ---
#alpha
png('figs/qq_alpha_gemma4.png', width=3.5, height=4, units='in', res=1200)
qq(a2$p_lrt,xaxt='n',yaxt='n',main='Chao1')
dev.off()

#beta1
png('figs/qq_beta1_gemma4.png', width=3.5, height=4, units='in', res=1200)
qq(b12$p_lrt,xaxt='n',yaxt='n',main='UniFrac PCoA axis 1')
dev.off()

#beta2
png('figs/qq_beta2_gemma4.png', width=3.5, height=4, units='in', res=1200)
qq(b22$p_lrt,xaxt='n',yaxt='n',main='UniFrac PCoA axis 2')
dev.off()

#beij
png('figs/qq_beij_gemma4.png', width=3.5, height=4, units='in', res=1200)
qq(beij$p_lrt,xaxt='n',yaxt='n',main='Beij pres/abs')
dev.off()

#micro
png('figs/qq_micro_gemma4.png', width=3.5, height=4, units='in', res=1200)
qq(micro$p_lrt,xaxt='n',yaxt='n',main='Micro pres/abs')
dev.off()

save.image('man_plots_gemma4.RData')


