setwd("~/Documents/projects/YRWA_GWAS_microbiome/angsd_fst")
library(qqman)
library(gtools)

load('plot_fst.RData')

#read in fst results
f<-read.table('sliding_window_fst.txt', skip=1)
colnames(f)<-c('region','chr','midPos','Nsites','fst')

summary(f$fst)
hist(f$fst,breaks=100)

#order df by chromosomes in alphanumeric order
f2<-f[mixedorder(f$chr),]
row.names(f2)<-1:nrow(f2)

#save chromosome names for plot
cnames<-unique(f2$chr)
#convert chromosomes to integer for qqman package
f2$chr<-as.integer(factor(f2$chr, levels = cnames)) #LEVELS have to be in right order, too
f2$snp<-rep('.',length(f2$chr)) #have to specify snp column for manhattan function!

#cut off scaffolds and mito (keep everything up to the Z=31)
f3<-f2[which(f2$chr %in% seq(1,31,by=1)),]
cnames3<-cnames[1:31]
png('man_fst_slidingwindow.png', width=6.5, height=2.5, units='in', res=1200)
manhattan(f3, chr = "chr", bp = "midPos", p = "fst", snp = "snp", chrlabs=substring(cnames3,4), ylim=c(0,1),
          ylab=expression(F[ST]),suggestiveline=F, logp=F,
          cex=0.5, cex.axis=0.6)
dev.off()

save.image('plot_fst.RData')
