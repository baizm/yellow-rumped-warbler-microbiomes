setwd("~/Documents/projects/YRWA_GWAS_microbiome/SNPRelate/")
#BiocManager::install("SNPRelate")

library(gdsfmt)
library(SNPRelate)
library(ggplot2)

load('make_pca_unique.RData')

#read-in the filtered vcf file
#this one was generated after re-making the joint genotypes with replicates removed
vcf<-('~/Documents/projects/YRWA_GWAS_microbiome/GBS_pipeline/exclude_replicates/biallelic_DP10_maf_thin_hwe_50.recode.vcf')

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
#change replicated id to match pca ID
metadata$id2[which(metadata$id2=='VF03T03.R1')]<-'VF03T03_1.1.R1'

#remove samples not in pca
x<-metadata$id2[(which(!(metadata$id2 %in% pca$sample.id)))] 
metadata3<-metadata[which(!(metadata$id2 %in% x)),] 
#arrange by id2
metadata3<-metadata3[order(metadata3$id2),]
#make sure they are the same:
identical(metadata3$id2,pca$sample.id) #TRUE!


ab <- data.frame(sample.id = pca$sample.id,
                 pop = metadata3$State,
                 proj = metadata3$Project,
                 EV1 = pca$eigenvect[,1],    # the first eigenvector
                 EV2 = pca$eigenvect[,2],    # the second eigenvector
                 stringsAsFactors = FALSE)

#plot pca
levels(ab$pop)<-list('1'='AB', '2'='AK','3'='AZ',
                     '4'='BC','5'='MI', '6'='NY',
                     '7'='PA','8'='UT')

levels(ab$proj)<-list('1'='ADK', '2'='AK_BC_myrtle','3'='AUWA_SW',
                     '4'='Opportunistic_MI','5'='Opportunistic_PA', '6'='YRWA_hybrid_zone_2019',
                     '7'='YRWA_hybrid_zone_2022')

#by geo
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
png('pca_unique.png', width=5.5, height=3.5, units='in', res=2400)
par(mar=c(5.1, 4.1, 2.1, 7.1), xpd=TRUE)
plot(ab$EV1, ab$EV2,pch=21,bg=alpha(colour = colors3[factor(ab$geo)], alpha=0.9),
     xlab="PC 1 (1.72%)", ylab="PC 2 (0.58%)",cex=1.2)
legend("topright", inset=c(-0.4,0), legend=c('Alaska','British Columbia','HZ 2019','HZ 2022',
                                             'Northeast','Southwest'), 
       pch=c(21), pt.bg=alpha(colour = colors3, alpha=0.9),title="Locality", 
       pt.cex=1.2,cex=0.7)
dev.off()

  

#calculate IBD values
ibd<-snpgdsIBDMoM(genofile, kinship=T,verbose=TRUE)
heatmap(ibd$kinship, Rowv = NA, Colv = NA, labRow = ibd$sample.id, labCol = ibd$sample.id)


ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)


hist(ibd.coeff$kinship)

ibd.coeff[which(ibd.coeff$kinship >0.01),] #one AK pair 

#identity-by-state
ibs <- snpgdsIBS(genofile, autosome.only=F)
heatmap(ibs$ibs, Rowv = NA, Colv = NA)

#write covariates file for GEMMA
#lat, long, pc1, pc2
covars<-data.frame(EV1=ab$EV1, EV2=ab$EV2, row.names=substr(ab$sample.id,1,7))
identical(row.names(metadata3), row.names(covars)) #TRUE, copy lat and long
covars$lat<-metadata3$Latitude
covars$long<-metadata3$Longitude
covars$year<-metadata3$Year
write.table(covars, 'covars.txt', sep=' ', quote=F, row.names=T, col.names=T)


save.image('make_pca_unique.RData')

