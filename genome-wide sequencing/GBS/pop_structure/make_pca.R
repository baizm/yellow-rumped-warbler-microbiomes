setwd("~/Documents/projects/YRWA_GWAS_microbiome/SNPRelate/")
#BiocManager::install("SNPRelate")

library(gdsfmt)
library(SNPRelate)

load('make_pca.RData')

#read-in the filtered vcf file
vcf<-('~/Documents/projects/YRWA_GWAS_microbiome/SNPRelate/biallelic_DP10_maf_thin_hwe_50.recode.vcf')

#reformat to gds and keep only bialleleic snps (no multiallelic or indels or structral vars)
snpgdsVCF2GDS(vcf,'gds',method='copy.num.of.ref') #use all sites coded as 0,1,2

snpgdsSummary('gds')

genofile <- snpgdsOpen('gds',allow.duplicate=T)

pca<-snpgdsPCA(genofile,autosome.only=F)

#calculate percent variation explained by first PCs
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#read in metadata file for locality information
metadata<-read.csv('~/Documents/projects/YRWA_GWAS_microbiome/metadata/metadata_yrwa_16s.tsv', sep='\t', header=T)
rownames(metadata)<-metadata$id
metadata<-metadata[which(metadata$Type=='Sample'),] #get rid of negatives
metadata$id2<-paste(metadata$id,"R1",sep='.') #make new column to match id in pca dataframe
#add the replicate sample
rep_ids<-pca$sample.id[grep('^VF03T03',pca$sample.id)] #pull out the id's of the replicate sample
metadata2<-rbind(metadata, metadata[rep('VF03T03', length(rep_ids)), ]) #add to bottom
metadata2$id2[(nrow(metadata2)-11):nrow(metadata2)]<-rep_ids #change the id to match pca
#remove samples not in pca
x<-metadata2$id2[(which(!(metadata2$id2 %in% pca$sample.id)))] 
metadata3<-metadata2[which(!(metadata2$id2 %in% x)),] 
#arrange by id2
metadata3<-metadata3[order(metadata3$id2),]
#make sure they are the same:
identical(metadata3$id2,pca$sample.id) #TRUE!


ab <- data.frame(sample.id = pca$sample.id,
                 pop = metadata3$State,
                 EV1 = pca$eigenvect[,1],    # the first eigenvector
                 EV2 = pca$eigenvect[,2],    # the second eigenvector
                 stringsAsFactors = FALSE)

#plot pca
levels(ab$pop)<-list('1'='AB', '2'='AK','3'='AZ','4'='BC','5'='MI', '6'='NY',
                     '7'='PA','8'='UT')
#colors<-c('gray','cornflowerblue','red','gray29','yellow','yellow4','orange','pink')
#plot(ab$EV1, ab$EV2, col=colors[factor(ab$pop)], xlab="PC 1 (3.51%)", ylab="PC 2 (1.68%)", 
#     pch=21,lwd=1.5,cex=0.9)
#legend("bottomleft", legend=c('AB','AK','AZ','BC','MI','NY','PA','UT'), 
 #      pch='o', col=colors, bty='n',cex=0.9)

#plot
ab$proj<-metadata3$Project
levels(ab$proj)<-list('1'='ADK', '2'='AK_BC_myrtle','3'='AUWA_SW',
                      '4'='Opportunistic_MI','5'='Opportunistic_PA', '6'='YRWA_hybrid_zone_2019',
                      '7'='YRWA_hybrid_zone_2022')

#create geography variable for plot
metadata3$geo<-rep('geo',length(metadata3$State))
metadata3$geo[which(metadata3$State=='NY' | metadata3$State=='MI' | metadata3$State=='PA')]<-'Northeast'
metadata3$geo[which(metadata3$State=='AK')]<-'Alaska'
metadata3$geo[which(metadata3$State=='BC' & metadata3$Project=='AK_BC_myrtle')]<-'British Columbia'
metadata3$geo[which(metadata3$State=='AZ' | metadata3$State=='UT')]<-'Southwest'
metadata3$geo[which(metadata3$Project=='YRWA_hybrid_zone_2019')]<-'HZ 2019'
metadata3$geo[which(metadata3$Project=='YRWA_hybrid_zone_2022')]<-'HZ 2022'

ab$geo<-metadata3$geo
levels(ab$geo)<-list('1'='Alaska','2'='British Columbia',
                      '3'='HZ 2019', '4'='HZ 2022', '5'='Northeast', '6'='Southwest')
colors3<-c('blue','cornflowerblue','gray29','gray','blue4','firebrick')
                          
#by geo
png('pca_replicates.png', width=5.5, height=3.5, units='in', res=2400)
par(mar=c(5.1, 4.1, 2.1, 7.1), xpd=TRUE)
plot(ab$EV1, ab$EV2,pch=21,bg=alpha(colour = colors3[factor(ab$geo)], alpha=0.9),
     xlab="PC 1 (3.51%)", ylab="PC 2 (1.68%)",cex=1.2)
legend("topright", inset=c(-0.4,0), legend=c('Alaska','British Columbia','HZ 2019','HZ 2022',
                                             'Northeast','Southwest'), 
       pch=c(21), pt.bg=alpha(colour = colors3, alpha=0.9),title="Locality", 
       pt.cex=1.2,cex=0.7)
dev.off()


#which are the outliers?
pca$sample.id[which(pca$eigenvect[,1] < -0.05)] #the replicates

#calculate IBD values
ibd<-snpgdsIBDMoM(genofile, kinship=T,verbose=TRUE)
ibd$kinship[151:165,151:165] #the section that contains replicates
heatmap(ibd$kinship, Rowv = NA, Colv = NA, labRow = ibd$sample.id, labCol = ibd$sample.id)

library(RColorBrewer)
png('ibd_replicates.png', width=5.5, height=3.5, units='in', res=2400)
heatmap(ibd$kinship,labRow = ibd$sample.id, labCol = ibd$sample.id,col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
legend(x="topright", legend=c("0","0.5"), 
       fill=colorRampPalette(brewer.pal(2, "Oranges"))(3))
dev.off()

heatmap(ibd$kinship[151:165,151:165], Rowv = NA, Colv = NA)

ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

hist(ibd.coeff$kinship)
summary(ibd.coeff$kinship) #mean=0.0009594

ibd.coeff[which(ibd.coeff$kinship >0.1 & ibd.coeff$kinship<0.5),] #one AK pair + replicate sample
ibd.vals.rep<-ibd.coeff[which(ibd.coeff$kinship >0.2 & ibd.coeff$kinship<0.5),] #only replicate samples
summary(ibd.vals.rep$kinship)


#identity-by-state
ibs <- snpgdsIBS(genofile, autosome.only=F)
heatmap(ibs$ibs, Rowv = NA, Colv = NA)

heatmap(ibs$ibs[151:165,151:165], Rowv = NA, Colv = NA)

save.image('make_pca.RData')

