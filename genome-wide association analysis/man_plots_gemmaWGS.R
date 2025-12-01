###---script for plotting results of GEMMA analysis of WGS data---
setwd("~/Documents/projects/YRWA_GWAS_microbiome/gemma/gemmaWGS/output/post_perm/")

library(qqman)
library(qvalue)
library(dplyr)
library(ggplot2)

load('man_plots_gemmaWGS.RData')

###-----alpha----
a<-read.table('../gemma_outWGS_a.assoc.txt.gz', header=T)

qa <- qvalue(p = a$p_lrt)
a$qvalue <- qa$qvalues
table(a$chr) 
unique(a$chr)

#re-order by chromosomes how I want -- already in order!, but still have to change to factor and make numeric...
chrs_ord<-unique(a$chr)
a <- a %>% arrange(factor(chr, levels = chrs_ord))
#copy and manipulate further
a2<-mutate(a, chr2 = as.numeric(factor(chr, levels = chrs_ord)))
#make numeric for plot
chrlabs <- c(1,"1a",2:4,'4a',5:15,17:29,"Z") #for plot

#make vector of sig snps to highlight
a2$rs2<-paste(rep('rs',length(a2$rs)),1:length(a2$rs),sep='') 
h_alpha<-a2$rs2[which(a2$p_lrt<1.672008e-08)] #threshold based on permutations
#plot!

png('figs/man_alpha_gemmaWGS.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(a2, chr = "chr2", bp = "ps", p = "p_lrt", snp = "rs2", chrlabs = chrlabs, ylim=c(0,8.5),
          highlight=h_alpha, main='Chao 1',suggestiveline=F, genomewideline = F,
          cex=0.5, cex.axis=0.6)
dev.off()

a2[which(a2$rs2 %in% h_alpha),] 

##----beta PCOA1-----
b1<-read.table('../gemma_outWGS_b1.assoc.txt.gz', header=T)

qb1 <- qvalue(p = b1$p_lrt)
b1$qvalue <- qb1$qvalues
table(b1$chr) 
unique(b1$chr)

#re-order by chromosomes how I want -- already in order!, but still have to change to factor and make numeric...
b1 <- b1 %>% arrange(factor(chr, levels = chrs_ord))
#copy and manipulate further
b12<-mutate(b1, chr2 = as.numeric(factor(chr, levels = chrs_ord))) #make numeric for plot

#make vector of sig snps to highlight
b12$rs2<-paste(rep('rs',length(b12$rs)),1:length(b12$rs),sep='') 
h_beta1<-b12$rs2[which(b12$p_lrt<3.075372e-08)] #threshold based on permutations

png('figs/man_beta1_gemmaWGS.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(b12, chr = "chr2", bp = "ps", p = "p_lrt", snp = "rs2", chrlabs = chrlabs, ylim=c(0,8.4),
          highlight=h_beta1, main='UniFrac PCoA 1',suggestiveline=F, genomewideline = F,
          cex=0.5, cex.axis=0.6)
dev.off()

b12[which(b12$rs2 %in% h_beta1),] 

##----beta PCOA2-----
b2<-read.table('../gemma_outWGS_b2.assoc.txt.gz', header=T)

qb2 <- qvalue(p = b2$p_lrt)
b2$qvalue <- qb2$qvalues
table(b2$chr) 
unique(b2$chr)

#re-order by chromosomes how I want -- already in order!, but still have to change to factor and make numeric...
b2 <- b2 %>% arrange(factor(chr, levels = chrs_ord))
#copy and manipulate further
b22<-mutate(b2, chr2 = as.numeric(factor(chr, levels = chrs_ord))) #make numeric for plot

#make vector of sig snps to highlight
b22$rs2<-paste(rep('rs',length(b22$rs)),1:length(b22$rs),sep='') 
h_beta2<-b22$rs2[which(b22$p_lrt<2.100484e-08)] #threshold based on permutations

png('figs/man_beta2_gemmaWGS.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(b22, chr = "chr2", bp = "ps", p = "p_lrt", snp = "rs2", chrlabs = chrlabs, ylim=c(0,9.8),
          highlight=h_beta2, main='UniFrac PCoA 2',suggestiveline=F, genomewideline = F,
          cex=0.5, cex.axis=0.6)
dev.off()

b22[which(b22$rs2 %in% h_beta2),] 

#---any overlaps?
identical(a2$rs2,b12$rs2); identical(b12$rs2,b22$rs2)
intersect(h_alpha,h_beta1)
intersect(h_beta1,h_beta2)
#no overlaps

#---- make QQ-plots ----
#alpha
png('figs/qq_alpha_wgs.png', width=3.5, height=4, units='in', res=1200)
qq(a2$p_lrt,xaxt='n',yaxt='n',main='Chao1')
dev.off()

#beta1
png('figs/qq_beta1_wgs.png', width=3.5, height=4, units='in', res=1200)
qq(b12$p_lrt,xaxt='n',yaxt='n',main='UniFrac PCoA axis 1')
dev.off()

#beta2
png('figs/qq_beta2_wgs.png', width=3.5, height=4, units='in', res=1200)
qq(b22$p_lrt,xaxt='n',yaxt='n',main='UniFrac PCoA axis 2')
dev.off()




save.image('man_plots_gemmaWGS.RData')
