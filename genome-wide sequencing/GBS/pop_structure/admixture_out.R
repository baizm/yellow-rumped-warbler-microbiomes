setwd("~/Documents/projects/YRWA_GWAS_microbiome/admixture_unique/")
library(dplyr)
load('admixture_out.RData')

#--------cross-validation error
cv<-read.table('outfiles/cv_k1-6.txt', sep=' ')

#CV plot for supplement
png('ev_error.png', width=3, height=3.5, units='in', res=2400)
plot(cv$V4, pch=21, col='black', bg='darkgray', cex=1.5,
     xlab='K', ylab='Cross-validation error')
lines(cv$V4, lty=2)
dev.off()

#------read in K=2 Qvalues
q<-read.table('outfiles/in_file_admixture.2.Q')

#get sample IDs from fam file
ids<-read.table('in_file_admixture.fam')
ids$V2<-substr(ids$V2, 1, 7)

#build data.frame to plot
rownames(q)<-ids$V2 #add ids

metadata<-read.csv('~/Documents/projects/YRWA_GWAS_microbiome/metadata/metadata_yrwa_16s.tsv', sep='\t', header=T)
rownames(metadata)<-metadata$id
metadata<-metadata[which(metadata$Type=='Sample'),] #get rid of negatives
#remove samples not in admixture file
x<-metadata$id[(which(!(metadata$id %in% rownames(q))))] 
metadata2<-metadata[which(!(metadata$id %in% x)),] 
#arrange by id
metadata2<-metadata2[order(metadata2$id),]
#make sure they are the same:
identical(metadata2$id,rownames(q)) #TRUE!

#create geography variable for plot
metadata2$geo<-rep('geo',length(metadata2$State))
metadata2$geo[which(metadata2$State=='NY' | metadata2$State=='MI' | metadata2$State=='PA')]<-'Northeast'
metadata2$geo[which(metadata2$State=='AK')]<-'Alaska'
metadata2$geo[which(metadata2$State=='BC' & metadata2$Project=='AK_BC_myrtle')]<-'British Columbia'
metadata2$geo[which(metadata2$State=='AZ' | metadata2$State=='UT')]<-'Southwest'
metadata2$geo[which(metadata2$Project=='YRWA_hybrid_zone_2019')]<-'HZ 2019'
metadata2$geo[which(metadata2$Project=='YRWA_hybrid_zone_2022')]<-'HZ 2022'
q$geo<-metadata2$geo
boxplot(q$V1~q$geo) #high=audubon, low=myrtle

colnames(q)<-c("Audubon's","Myrtle",'geo')

#order q by populations and ancestry
q$geo<-factor(q$geo, levels=c('Southwest','HZ 2019','HZ 2022','British Columbia','Alaska','Northeast'))

boxplot(q$`Audubon's`~q$geo, outline=F, xlab='', ylab='Q')
stripchart(q$`Audubon's`[which(q$geo=='Southwest')],vertical=T,pch=21,bg='firebrick', method='jitter',add=T)
stripchart(q$`Audubon's`[which(q$geo=='HZ 2019')],vertical=T,pch=21,bg='gray29', method='jitter',add=T, at=2)
stripchart(q$`Audubon's`[which(q$geo=='HZ 2022')],vertical=T,pch=21,bg='gray', method='jitter',add=T, at=3)
stripchart(q$`Audubon's`[which(q$geo=='British Columbia')],vertical=T,pch=21,bg='cornflowerblue', method='jitter',add=T, at=4)
stripchart(q$`Audubon's`[which(q$geo=='Alaska')],vertical=T,pch=21,bg='blue', method='jitter',add=T, at=5)
stripchart(q$`Audubon's`[which(q$geo=='Northeast')],vertical=T,pch=21,bg='blue4', method='jitter',add=T, at=6)

#arrange by geo, then by Q
q2<-q %>% 
    arrange(factor(geo, levels = c('Southwest','HZ 2019','HZ 2022','British Columbia','Alaska','Northeast')),
            desc(`Audubon's`))

table(q2$geo)

#structure plot
png('barplot_k2.png', width=6, height=3, units='in', res=2400)
bp<-barplot(t(as.matrix(q2[,1:2])), col=c('firebrick','dodgerblue3'), ylim=c(0,1),
        xlab="", ylab="Admixture proportions", border=NA, space=0, cex.axis=1)
abline(v=bp[c(33, 33+51, 33+51+67, 33+51+67+34, 33+51+67+34+40)], lwd=1.5)
dev.off()

#----test for substructure w/in Myrtle?
allo_ids<-rownames(q)[which(q$geo %in% c('Northeast','Alaska','British Columbia'))]
#add .R1 to match ids in VCF file
allo_ids<-paste(allo_ids,'.R1', sep='')
write.table(allo_ids, 'substructure_myrtles/allo_ids.txt', quote=F, sep='\n', col.names=F, row.names=F) #write to file for VCF filtering
#ran ADMIXTURE on CCR, return here to look at results!

#cross-validation error
cv_myr<-read.table('substructure_myrtles/outfiles/cv_k1-4.txt', sep=' ')
#CV plot for supplement?
png('substructure_myrtles/ev_error_myr.png', width=3, height=3.5, units='in', res=2400)
plot(cv_myr$V4, pch=21, col='black', bg='darkgray', cex=1.5,
     xlab='K', ylab='Cross-validation error')
lines(cv_myr$V4, lty=2)
dev.off()


#read in K=3 Qvalues
qm<-read.table('substructure_myrtles/outfiles/myrtles_admixture.3.Q')

#get sample IDs from fam file
idsm<-read.table('substructure_myrtles/myrtles_admixture.fam')
idsm$V2<-substr(idsm$V2, 1, 7)

#build data.frame to plot
rownames(qm)<-idsm$V2 #add ids

#grab geo column in metadata2 data.frame
test<-metadata2[which(metadata2$id %in% rownames(qm)),]
identical(test$id, rownames(qm))
qm$geo<-test$geo

qm2<-qm %>% 
  arrange(factor(geo, levels = c('British Columbia','Alaska','Northeast')),
          desc(V1))

table(qm2$geo)
str(qm2$geo)

#structure plot
png('substructure_myrtles/barplot_k3_myr.png', width=6, height=3.1, units='in', res=2400)
bpm<-barplot(t(as.matrix(qm2[,1:3])), col=c('cornflowerblue','blue','blue4'), ylim=c(0,1),
            xlab="", ylab="Admixture proportions", border=NA, space=0., cex.axis = 1)
abline(v=bpm[c(34.5, 34.5+40)], lwd=1.5,col='black')
dev.off()

save.image('admixture_out.RData')

#write Q dataframe to file for microbiome script
saveRDS(q, 'q.RDS')

#write IDs of myr and aud for fst analyses
myr_ids<-rownames(q[which(q$`Audubon's` <0.05),])
myr_ids<-paste('/projects/academic/mbaiz/YRWA_GBS_2023/bams_sorted2/',myr_ids,'.R1.sorted.bam', sep='')
write.table(myr_ids,'myr_ids.txt', quote=F,sep='/n',row.names=F,col.names=F)

aud_ids<-rownames(q[which(q$`Audubon's` >0.95),])
aud_ids<-paste('/projects/academic/mbaiz/YRWA_GBS_2023/bams_sorted2/',aud_ids,'.R1.sorted.bam', sep='')
write.table(aud_ids,'aud_ids.txt', quote=F,sep='/n',row.names=F,col.names=F)
