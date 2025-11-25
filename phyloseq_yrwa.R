setwd("~/Documents/projects/YRWA_GWAS_microbiome/qiime/phyloseq/")
library(phyloseq)
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(vegan)
library(reshape2)
library(decontam)
library(gridExtra)
library(plyr)
library(lubridate)
library(maps)
library(mapdata)
library(maptools)
data(wrld_simpl)
library(btools) #for faiths pd
library(Maaslin2)
library(geodist)
library(car)

load('phyloseq_yrwa.RData')

#read in post-decontam filtered data from qiime2
########-------create physeq object---------------------########
####------read in OTU table, has to be a matrix---------------------
otu_table<-read.csv('../qiime_out/feature-table.tsv', sep='\t', header=T, skip=1, check.names = F)
#convert to matrix
otumat<-as.matrix(otu_table[,2:ncol(otu_table)])
#add rownames that are OTU id
rownames(otumat)<-otu_table[,1]
###-------read in taxonomy table-----------------
tax_table<-read.csv('../qiime_out/decontam/silva_taxonomy.tsv', sep='\t', header=F, skip=2)
tax_table2<-separate(data = tax_table, col = V2, 
                     into = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = "; ")
colnames(tax_table2)[1]<-'Feature.ID'
colnames(tax_table2)[9]<-'Confidence'
#it wasn't filtered, so i'll have to delete mt,cp,un,euk,decontam...
tax_table2<-tax_table2[which(tax_table2$Feature.ID %in% otu_table$`#OTU ID`),] # OTUs after filtering all mt,cp,un,euk,decontam
#convert to matrix
taxmat<-as.matrix(tax_table2[,2:8])
rownames(taxmat)<-tax_table2$Feature.ID
###-------read in sample data---------------------
sampledata<-read.csv('../../metadata/metadata_yrwa_16s.tsv', sep='\t', header=T,)
rownames(sampledata)<-sampledata$id
x<-sampledata$id[(which(!(sampledata$id %in% colnames(otumat))))]
sampledata<-sampledata[which(!(sampledata$id %in% x)),] #get rid of samples filtered by qiime (all negatives)
sampledata<-sampledata[,2:ncol(sampledata)]
###-------read in nwk tree with ape package-----------
tree<-read.tree('../qiime_out/tree.nwk')
###----------combine into phyloseq object-------------
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(sampledata)
physeq = phyloseq(OTU, TAX, SAM, tree)
#sums of OTUs across samples

##----housekeeping------
table(physeq@sam_data$Type, useNA = 'always') 
#pull out all but negatives
exl_negs<-rownames(physeq@sam_data)[which(!(grepl('Neg', sample_data(physeq)$Type)))]
physeq2<-prune_samples(rownames(physeq@sam_data) %in% exl_negs, physeq) 
physeq2<-prune_taxa(taxa_sums(physeq2) > 0, physeq2) #get rid of zero-sum OTUs
table(physeq2@sam_data$Species, useNA='always') 
physeq2 

#id recaps
table(duplicated(physeq2@sam_data$Band.)) #1 duplicated

#write sample names to file for GBS pipeline (genotyping)
test<-paste(rownames(physeq2@sam_data), "*", sep='')
test2<-paste(paste(as.character(test), sep="' '"), collapse=' ') #use this one
write.table(test2, 'sample_ids_genotyping.txt',quote=F, row.names=F,col.names=F)
#in gVCF folder, run
#cp $(<sample_ids_genotyping.txt) ../gVCFs2

####---------find rarefaction depth-------------
#make data frame with total read counts after filtering using sample_sums
sdt = data.table(as(sample_data(physeq2), "data.frame"),
                 TotalReads = sample_sums(physeq2), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
summary(sdt$TotalReads) #range 
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
#look at rarecurve
tab <- otu_table(physeq2)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
pdf(file='rarecurve.pdf', height=4, width=4)
rarecurve(tab, step=50, lwd=0.8, label=F, xlab='Read count',
          ylab='ASVs',xlim=c(0,7000))
dev.off()

#how many samples below 3000 reads?
dim(sdt[which(sdt$TotalReads < 3000),]) #78 -- too many
dim(sdt[which(sdt$TotalReads < 2000),]) #64 
dim(sdt[which(sdt$TotalReads < 1600),]) #55 
dim(sdt[which(sdt$TotalReads < 1500),]) #52 
dim(sdt[which(sdt$TotalReads < 1000),]) #44 

table(sdt$Project[which(sdt$TotalReads>2000)],
      sdt$Year[which(sdt$TotalReads>2000)])

table(sdt$Project[which(sdt$TotalReads>1600)],
      sdt$Year[which(sdt$TotalReads>1600)]) #-- can afford to be a little more conservative

ps<-rarefy_even_depth(physeq2,sample.size=2036,rngseed = 999, replace = F, trimOTUs = T, verbose = T) #threshold is below 2000
ps@sam_data$Year<-as.factor(ps@sam_data$Year)

#------remove recap-----####
ps<-prune_samples(!(rownames(ps@sam_data) %in% c('TF18T11')), ps) 
ps<-prune_taxa(taxa_sums(ps) > 0, ps) #get rid of zero-sum OTUs
ps 
table(duplicated(ps@sam_data$Band.)) #check!
table(ps@sam_data$Project, ps@sam_data$Species)

###----map loaclities----###
summary(ps@sam_data$Latitude)
summary(ps@sam_data$Longitude)
#read in .shp files downloaded from IUCN, and plot it
au_shp<-readShapePoly('../../map/AUWA_redlist_species_data_70b19554-c38e-4188-bc9f-1328b0df842a/data_0.shp')
my_shp<-readShapePoly('../../map/MYWA_redlist_species_data_c6ddfdf5-51f8-4213-af83-f1fcd4005650/data_0.shp')

png('figs/sampling_map_16S_v2.png', height=3.5, width=4, units='in', res=2400)
#png('figs/sampling_map_16S_nopoints.png', height=3.5, width=4, units='in', res=2400)
map('worldHires',c('USA', 'Canada','Mexico'), col='gray93', fill=T, lwd=0.1, xlim=c(-170,-50),ylim=c(25,71))
map.axes(cex.axis=0.7)
map.scale(cex=0.4, ratio=F, x=-165, y=35)
#add ranges
str(au_shp$SEASONAL)
palette(c('firebrick','firebrick','gray93','gray93','gray93')) #resident, breeding, .., .., ..
plot(au_shp, add=T, col=factor(au_shp$SEASONAL), border=F)
#plot(au_shp, add=T, col=alpha('firebrick', 0.8), border=F)
str(my_shp$SEASONAL)
palette(c('dodgerblue3','gray93','gray93','gray93','gray93')) #breeding, non-breeding,...,...,...
plot(my_shp, add=T, col=factor(my_shp$SEASONAL), border=F) 
#plot(my_shp, add=T, col=alpha('dodgerblue3', 0.6), border=F)
#all samples
points(ps@sam_data$Longitude, ps@sam_data$Latitude, pch=21, bg='black',col='white', lwd=.6, cex=0.6)
dev.off()

points(ps@sam_data$Longitude, ps@sam_data$Latitude, pch=21, bg=factor(ps@sam_data$Year), lwd=.4)
points(ps@sam_data$Longitude, ps@sam_data$Latitude, pch=21, bg=factor(ps@sam_data$Project), lwd=.4)

#----read in plumage scores and other extra data-------
p1<-read.csv('../../metadata/w:plumage/2022_YRWA_banding_records_ab.xlsx - Sheet1.csv', header=T)
rownames(p1)<-p1$Record.Number
p1<-p1[,c(32:40,28,31,30)]
p2<-read.csv('../../metadata/w:plumage/yrwa_data_2019_re-ordered.csv', header=T)
rownames(p2)<-p2$ID
p2<-p2[,c(19:27,34:36)]
p3<-rbind(p1,p2) #combine
p3<-p3[order(rownames(p3)),] #order
p3<-p3[which(rownames(p3) %in% rownames(ps@sam_data)),] #keep only birds in ps, N=124
test<-merge(ps@sam_data, p3, by='row.names', all=T) #merge
#bill_exposed
b<-read.csv('../../metadata/w:plumage/2022_YRWA_banding_records_ab.xlsx - Sheet1.csv',header=T)
b<-data.frame(Bill_Exposed=b$Bill_Exposed, id=b$Record.Number)
b<-b[which(b$id %in% rownames(ps@sam_data)),]
test2<-merge(test,b, by.x='Row.names', by.y='id', all=T) #merge
#weight
w<-read.csv('../../metadata/w:plumage/yrwa_data_2019_re-ordered.csv', header=T)
rownames(w)<-w$ID
w<-w[,c(37:39)]
w<-w[order(rownames(w)),]
w<-w[which(rownames(w) %in% rownames(ps@sam_data)),]
test3<-merge(test2,w, by.x='Row.names', by.y='row.names', all=T) #final merge
dim(test3); dim(ps@sam_data) #check
identical(as.character(rownames(ps@sam_data)), as.character(test3$Row.names)) #TRUE
#replace with new data
new_sampledata<-cbind(ps@sam_data, test3[c(29:44)]) #write new sampledata
sample_data(ps)<-new_sampledata
#fix certain column types
ps@sam_data$extraction.year<-factor(ps@sam_data$extraction.year)
ps@sam_data$Wing_patch<-as.integer(ps@sam_data$Wing_patch)

hist(ps@sam_data$Auricular_extent)
hist(ps@sam_data$Throat_extent_yellow, breaks=20)
hist(ps@sam_data$Throat_saturation_yellow, breaks=20)
hist(ps@sam_data$Wing_patch)

#calculate plumage score (high scores=MYWA, low scores=AUWA)
ps@sam_data$plum_score<-((100-ps@sam_data$Throat_extent_yellow)+(100-ps@sam_data$Throat_saturation_yellow)+
        (ps@sam_data$Throat_white_corners)+(ps@sam_data$Auricular_extent)+
        (ps@sam_data$Auricular_saturation)+(ps@sam_data$SPOT)+
        (ps@sam_data$LINE)+(100-ps@sam_data$Wing_patch))/8
hist(ps@sam_data$plum_score, breaks=50)
hist(ps@sam_data$plum_score[which(ps@sam_data$Project=='YRWA_hybrid_zone_2022')])
hist(ps@sam_data$plum_score[which(ps@sam_data$Project=='YRWA_hybrid_zone_2019')])


#------calculate alpha values----------
alpha_r<-estimate_richness(ps, measures=c('Shannon','Observed','Chao1'))
alpha_r$PD<-estimate_pd(ps)[,1] #grab first column that has faith's pd (second column is species richness--#otus)
alpha_r$Species<-ps@sam_data$Species
alpha_r$project<-ps@sam_data$Project
alpha_r$elevation<-ps@sam_data$Elevation
alpha_r$sex<-ps@sam_data$Sex
alpha_r$age<-ps@sam_data$Age
alpha_r$Year<-ps@sam_data$Year
alpha_r$extraction.year<-ps@sam_data$extraction.year
alpha_r$plum_score<-ps@sam_data$plum_score
alpha_r$latitude<-ps@sam_data$Latitude
alpha_r$longitude<-ps@sam_data$Longitude
alpha_r$state<-ps@sam_data$State


#--------analyses with indiviuals that have both admixture Q and 16S data ----
q<-readRDS('~/Documents/projects/YRWA_GWAS_microbiome/admixture_unique/q.RDS')

#subset to include only birds with Q
psq<-prune_samples(rownames(ps@sam_data)[which(rownames(ps@sam_data) %in% rownames(q))], ps) #microbiome birds with Q data
psq<-prune_taxa(taxa_sums(psq) > 0, psq) #get rid of zero-sum OTUs

qq<-q[which(rownames(q) %in% rownames(psq@sam_data)),] #Q birds w/microbiome data

identical(rownames(psq@sam_data),rownames(qq)) #TRUE

psq@sam_data$aud<-qq$`Audubon's`
psq@sam_data$geo<-qq$geo
psq 
table(psq@sam_data$geo, psq@sam_data$Year)
table(psq@sam_data$geo)
table(psq@sam_data$Year, psq@sam_data$qsp)

#assign Qspecies
psq@sam_data$qsp<-rep('qsp',length(psq@sam_data$aud)) #make a new column, then add variables
psq@sam_data$qsp[which(psq@sam_data$aud > 0.95)]<-'aud'
psq@sam_data$qsp[which(psq@sam_data$aud < 0.05)]<-'myr'
psq@sam_data$qsp[which(psq@sam_data$aud < 0.95 & psq@sam_data$aud > 0.05 )]<-'hyb'
table(psq@sam_data$qsp, useNA = 'always')

#lump hz transects for geographic variable
geo2<-as.vector(psq@sam_data$geo) #copy
geo2[which(geo2=='HZ 2019' | geo2=='HZ 2022')]<-'Hybrid Zone'
psq@sam_data$geo2<-geo2

table(psq@sam_data$qsp, psq@sam_data$Year, useNA = 'always')

ordinate(psq_sex, "PCoA", "unifrac") %>%
  plot_ordination(psq_sex, ., color='qsp', title = "UniFrac", type='samples') +
  stat_ellipse(aes(color=qsp, group=qsp), level=0.5, lty=2, lwd=0.4) +
  geom_point(aes(fill=qsp, shape=qsp), size=2, color='black') +
  scale_color_manual(values=c('firebrick','darkgray','dodgerblue')) + 
  scale_fill_manual(values=alpha(c('firebrick','darkgray','dodgerblue'), 0.6)) +
  scale_shape_manual(values=c(21,21,21)) +
  theme(axis.text = element_text(size = 14)) +
  theme_bw() + labs(colour="Admixture group") 

png('figs/pcoa_psq_sex.png', width=4, height=3, units='in', res=2400)
ordinate(psq_sex, "PCoA", "bray") %>%
  plot_ordination(psq_sex, ., color='qsp', title = "Bray Curtis", type='samples') +
  stat_ellipse(aes(color=qsp, group=qsp), level=0.9, lty=2, lwd=0.4) +
  geom_point(aes(fill=qsp, shape=qsp), size=2, color='black') +
  scale_color_manual(values=c('firebrick','darkgray','dodgerblue')) + 
  scale_fill_manual(values=alpha(c('firebrick','darkgray','dodgerblue'), 0.6)) +
  scale_shape_manual(values=c(21,21,21)) +
  theme(axis.text = element_text(size = 14)) +
  theme_bw() + labs(colour="Admixture group") 
dev.off()



#----test variables for significance ----
table(psq@sam_data$Sex, useNA = 'always') 
#get rid of ?sex
psq@sam_data$Sex[which(psq@sam_data$Sex=='F?')]<-NA
psq@sam_data$Sex[which(psq@sam_data$Sex=='M?')]<-NA
table(psq@sam_data$Sex, psq@sam_data$Year) 

#make ps object with NA sex removed so we can test using adonis2
psq_sex<-prune_samples(rownames(psq@sam_data)[which(psq@sam_data$Sex %in% c('M','F'))], psq) 
psq_sex<-prune_taxa(taxa_sums(psq_sex) > 0, psq_sex) #get rid of zero-sum OTUs

bray_full<-phyloseq::distance(psq_sex, method='bray')
jac_full<-phyloseq::distance(psq_sex, method='jaccard', binary=T)
uni_full<-phyloseq::distance(psq_sex, method='unifrac')
wuni_full<-phyloseq::distance(psq_sex, method='wunifrac')

adonis2(bray_full ~ Year+qsp+Sex+geo2, data=data.frame(psq_sex@sam_data), by='margin') 
adonis2(jac_full ~ Year+qsp+Sex+geo2, data=data.frame(psq_sex@sam_data), by='margin') 
adonis2(uni_full ~ Year+qsp+Sex+geo2, data=data.frame(psq_sex@sam_data), by='margin') 
adonis2(wuni_full ~ Year+qsp+Sex+geo2, data=data.frame(psq_sex@sam_data), by='margin') 

permutest(betadisper(bray_full, psq_sex@sam_data$Year)) 
permutest(betadisper(bray_full, psq_sex@sam_data$qsp)) 
permutest(betadisper(bray_full, psq_sex@sam_data$Sex)) 
permutest(betadisper(bray_full, psq_sex@sam_data$geo2)) 
permutest(betadisper(jac_full, psq_sex@sam_data$Year)) 
permutest(betadisper(jac_full, psq_sex@sam_data$qsp)) 
permutest(betadisper(jac_full, psq_sex@sam_data$Sex)) 
permutest(betadisper(jac_full, psq_sex@sam_data$geo2)) 
permutest(betadisper(uni_full, psq_sex@sam_data$Year)) 
permutest(betadisper(uni_full, psq_sex@sam_data$qsp)) 
permutest(betadisper(uni_full, psq_sex@sam_data$Sex)) 
permutest(betadisper(uni_full, psq_sex@sam_data$geo2)) 
permutest(betadisper(wuni_full, psq_sex@sam_data$Year)) 
permutest(betadisper(wuni_full, psq_sex@sam_data$qsp)) 
permutest(betadisper(wuni_full, psq_sex@sam_data$Sex)) 
permutest(betadisper(wuni_full, psq_sex@sam_data$geo2)) 


alpha_sex<-alpha_r[which(rownames(alpha_r) %in% rownames(psq_sex@sam_data)),]
identical(rownames(psq_sex@sam_data),rownames(alpha_sex)) #TRUE
alpha_sex$aud<-psq_sex@sam_data$aud
alpha_sex$qsp<-psq_sex@sam_data$qsp
alpha_sex$geo2<-psq_sex@sam_data$geo2

lm_chao1_full <- lm(Chao1 ~ Year+qsp+sex+geo2, data=alpha_sex)
lm_chao1_full %>% car::Anova() 
lm_Shannon_full <- lm(Shannon ~ Year+qsp+sex+geo2, data=alpha_sex)
lm_Shannon_full %>% car::Anova() 
lm_pd_full <- lm(PD ~ Year+qsp+sex+geo2, data=alpha_sex)
lm_pd_full %>% car::Anova() 


#---write files for GEMMA
#log transform alpha for gemma
log_chao1<-data.frame(log_chao1=log(alpha_r$Chao1), row.names=row.names(alpha_r)) #includes sex unknown
log_chao1<-log_chao1[which(rownames(log_chao1) %in% rownames(psq@sam_data))]
write.table(log_chao1, 'log_chao1.txt', sep=' ', quote=F, row.names=T, col.names=T)

#----maaslin2
maas_all_qsp = Maaslin2(
  input_data = data.frame(psq_sex@otu_table), 
  input_metadata = data.frame(psq_sex@sam_data), 
  output = "maaslin2_output/psq_sex_all_qsp", 
  fixed_effects = c("Year","qsp","geo2","Sex"),
  reference = c("Year,2017","qsp,myr","geo2,Northeast"))

sig_asv_all<-data.frame(maas_all_qsp$results)
sig_asv_all<-sig_asv_all[which(sig_asv_all$qval < 0.25),]
sig_asv_all$feature<-sub('X','',sig_asv_all$feature) 
tax_for_sig_asv_all<-data.frame(psq_sex@tax_table)
tax_for_sig_asv_all<-tax_for_sig_asv_all[which(rownames(tax_for_sig_asv_all) %in% sig_asv_all$feature),]

sig_asv_all<-merge(sig_asv_all, tax_for_sig_asv_all, by.x='feature', by.y='row.names', all.x=T)
write.csv(sig_asv_all, 'maaslin2_results_w_taxonomy.csv')

length(unique(sig_asv_all$feature))
sig_asv_all[which(duplicated(sig_asv_all$feature)),]


#make plots of significant asvs
asvs_to_plot<-psq_sex@otu_table[which(rownames(psq_sex@otu_table) %in% sig_asv_all$feature),]
asvs_to_plot_counts<-data.frame(psq_sex@sam_data)
#two qsp asvs
asvs_to_plot_counts$count_453=as.vector(asvs_to_plot['453cf07d3e96e304f23fbde7b355fe70',])
asvs_to_plot_counts$count_947=as.vector(asvs_to_plot['9477d45e83a2f8d1f61fd3cb9db73a0f',])
asvs_to_plot_counts$count_90f=as.vector(asvs_to_plot['90f44084459949dd4fd965e1fc355026',])
asvs_to_plot_counts$count_bcb=as.vector(asvs_to_plot['bcb2bfaaf593f54d9fbef25dae1c88e8',])
asvs_to_plot_counts$count_b8c=as.vector(asvs_to_plot['b8c4ed03d7d91546b9f6ad80b80135a3',])


png(file='figs/counts_micro.png', height=4, width=3, units='in', res=2400)
boxplot(asvs_to_plot_counts$count_453~asvs_to_plot_counts$qsp, xlab='', ylab='Counts of Microbacteriaceae ASV',
        main='', boxwex=0.5, outline=F,ylim=c(0,140),col='white')
stripchart(asvs_to_plot_counts$count_453[which(asvs_to_plot_counts$qsp=='aud')],vertical=T,pch=21,bg=alpha('firebrick',0.7), method='jitter',add=T)
stripchart(asvs_to_plot_counts$count_453[which(asvs_to_plot_counts$qsp=='hyb')],vertical=T,pch=21,bg=alpha('gray83', 0.7), method='jitter',add=T,at=2)
stripchart(asvs_to_plot_counts$count_453[which(asvs_to_plot_counts$qsp=='myr')],vertical=T,pch=21,bg=alpha('dodgerblue', 0.7), method='jitter',add=T,at=3)
dev.off()

png(file='figs/counts_beij.png', height=4, width=3, units='in', res=2400)
boxplot(asvs_to_plot_counts$count_947~asvs_to_plot_counts$qsp, xlab='', ylab='Counts of Beijerinckiaceae ASV',
        main='', boxwex=0.5, outline=F,ylim=c(0,170),col='white')
stripchart(asvs_to_plot_counts$count_947[which(asvs_to_plot_counts$qsp=='aud')],vertical=T,pch=21,bg=alpha('firebrick',0.7), method='jitter',add=T)
stripchart(asvs_to_plot_counts$count_947[which(asvs_to_plot_counts$qsp=='hyb')],vertical=T,pch=21,bg=alpha('gray83', 0.7), method='jitter',add=T,at=2)
stripchart(asvs_to_plot_counts$count_947[which(asvs_to_plot_counts$qsp=='myr')],vertical=T,pch=21,bg=alpha('dodgerblue', 0.7), method='jitter',add=T,at=3)
dev.off()


#check counts 
psq_sex@otu_table[which(rownames(psq_sex@otu_table)=='bcb2bfaaf593f54d9fbef25dae1c88e8')] #not in SW, stating at UF15
psq_sex@otu_table[which(rownames(psq_sex@otu_table)=='90f44084459949dd4fd965e1fc355026')] #not in SW, stating at UF15

#geo2 associations
png(file='figs/counts_beij_geo2.png', height=4, width=3.5, units='in', res=2400)
boxplot(asvs_to_plot_counts$count_90f~asvs_to_plot_counts$geo2, xlab='', ylab='Counts of Beijerinckiaceae ASV',
        main='', boxwex=0.5, outline=F,ylim=c(0,80),col='white')
stripchart(asvs_to_plot_counts$count_90f[which(asvs_to_plot_counts$geo2=='Southwest')],vertical=T,pch=21,bg=alpha('black',0.4), method='jitter',add=T)
stripchart(asvs_to_plot_counts$count_90f[which(asvs_to_plot_counts$geo2=='Hybrid Zone')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=2)
stripchart(asvs_to_plot_counts$count_90f[which(asvs_to_plot_counts$geo2=='British Columbia')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=3)
stripchart(asvs_to_plot_counts$count_90f[which(asvs_to_plot_counts$geo2=='Alaska')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=4)
stripchart(asvs_to_plot_counts$count_90f[which(asvs_to_plot_counts$geo2=='Northeast')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=5)
dev.off()

png(file='figs/counts_unclass.png', height=4, width=3.5, units='in', res=2400)
boxplot(asvs_to_plot_counts$count_bcb~asvs_to_plot_counts$geo2, xlab='', ylab='Counts of unclassified ASV',
        main='', boxwex=0.5, outline=F,ylim=c(0,110),col='white')
stripchart(asvs_to_plot_counts$count_bcb[which(asvs_to_plot_counts$geo2=='Southwest')],vertical=T,pch=21,bg=alpha('black',0.4), method='jitter',add=T)
stripchart(asvs_to_plot_counts$count_bcb[which(asvs_to_plot_counts$geo2=='Hybrid Zone')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=2)
stripchart(asvs_to_plot_counts$count_bcb[which(asvs_to_plot_counts$geo2=='British Columbia')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=3)
stripchart(asvs_to_plot_counts$count_bcb[which(asvs_to_plot_counts$geo2=='Alaska')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=4)
stripchart(asvs_to_plot_counts$count_bcb[which(asvs_to_plot_counts$geo2=='Northeast')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=5)
dev.off()

png(file='figs/counts_blasto.png', height=4, width=3.5, units='in', res=2400)
boxplot(asvs_to_plot_counts$count_b8c~asvs_to_plot_counts$geo2, xlab='', ylab='Counts of Blastococcus ASV',
        main='', boxwex=0.5, outline=F,ylim=c(0,200),col='white')
stripchart(asvs_to_plot_counts$count_b8c[which(asvs_to_plot_counts$geo2=='Southwest')],vertical=T,pch=21,bg=alpha('black',0.4), method='jitter',add=T)
stripchart(asvs_to_plot_counts$count_b8c[which(asvs_to_plot_counts$geo2=='Hybrid Zone')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=2)
stripchart(asvs_to_plot_counts$count_b8c[which(asvs_to_plot_counts$geo2=='British Columbia')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=3)
stripchart(asvs_to_plot_counts$count_b8c[which(asvs_to_plot_counts$geo2=='Alaska')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=4)
stripchart(asvs_to_plot_counts$count_b8c[which(asvs_to_plot_counts$geo2=='Northeast')],vertical=T,pch=21,bg=alpha('black', 0.4), method='jitter',add=T,at=5)
dev.off()


#--pcoa eigens by hybrid index
#pull out PCoA values for gemma
ord2<-ordinate(psq_sex, "PCoA", "unifrac")
ord2$values$Relative_eig[1:2] #0.08609996 0.06861285 (amount of variation explained)
pcoa_ax2<-data.frame(pcoa1=ord2$vectors[,1], pcoa2=ord2$vectors[,2], row.names=row.names(ord2$vectors))
hist(pcoa_ax2$pcoa1); hist(pcoa_ax2$pcoa2) #pretty normal looking
#add q values to df
identical(rownames(psq_sex@sam_data), rownames(pcoa_ax2))
pcoa_ax2$aud<-psq_sex@sam_data$aud
cor.test(pcoa_ax2$pcoa2, pcoa_ax2$aud) #t = 6.302, df = 185, p-value = 2.096e-09, cor=0.4920397
cor.test(pcoa_ax2$pcoa1, pcoa_ax2$aud) #ns


##--- top taxa----- use ps n=222, including those w/o q values or sex known
topp<-as.data.frame(ps@tax_table)
table(topp$Phylum, useNA = 'always')
sort(table(topp$Phylum), decreasing=T)
sort(table(topp$Phylum)/sum(table(topp$Phylum)), decreasing=T) #top phyla by proportion

toppn<-names(sort(table(topp$Phylum), decreasing=T)[1:8]) #top 8 phyla

ps@sam_data$geo2<-qq$geo[match(rownames(ps@sam_data), rownames(qq))]
ps@sam_data$geo2[grep("^C", rownames(ps@sam_data))]<-'HZ 2019'
ps@sam_data$geo2[grep("^R", rownames(ps@sam_data))]<-'Northeast'
ps@sam_data$geo2[rownames(ps@sam_data)=='UF07T05']<-'Northeast'
ps@sam_data$geo2[rownames(ps@sam_data)=='VE31T02']<-'HZ 2022'
ps@sam_data$geo2[rownames(ps@sam_data) %in% c('VF09Z02','VF12Z04','VF13J01','VF14Z03','VF19J02','VF21Z03','VF25Z03')]<-'British Columbia'


glom<-tax_glom(ps, taxrank = "Phylum", NArm=F) #glom by phylum
mglom<-psmelt(glom)
mglom$Phylum[which(!(mglom$Phylum %in% toppn))]<-'other' #change non-common phyla to other 
table(mglom$Phylum)

png('figs/supp/relative_abundance_yrwa.png', width=10, height=7, units='in', res=1200)
ggplot(mglom, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill", width=1) +
  theme(axis.title.x=element_blank()) +
  theme(panel.background = element_blank()) +
  scale_fill_manual(values=c('gray','#F64C32','#FF8A25','#FFD75F','#0078FF','#09A4DB','#006E9A','#008D9A','#465D98'))+
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  facet_wrap(~geo2, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(strip.background=element_rect(colour="black",fill="white")) +
  theme(panel.spacing = unit(0.7, "lines"))
dev.off()


#map for supplement of WGS sample localities
png('figs/sampling_map_16S_v3.png', height=3.5, width=3, units='in', res=2400)
#map('worldHires',c('USA', 'Canada','Mexico'), col='gray93', fill=T, lwd=0.1, xlim=c(-170,-50),ylim=c(25,71))
map('worldHires',c('USA', 'Canada','Mexico'), col='gray93', fill=T, lwd=0.1, xlim=c(-123,-112),ylim=c(48,56.5))
map.axes(cex.axis=0.7)
#add ranges
str(au_shp$SEASONAL)
palette(c('firebrick','firebrick','gray93','gray93','gray93')) #resident, breeding, .., .., ..
plot(au_shp, add=T, col=factor(au_shp$SEASONAL), border=F)
str(my_shp$SEASONAL)
palette(c('dodgerblue3','gray93','gray93','gray93','gray93')) #breeding, non-breeding,...,...,...
plot(my_shp, add=T, col=factor(my_shp$SEASONAL), border=F) 
#only WGS samples versus non-wgs samples
wgs_ids<-read.table('~/Documents/projects/YRWA_GWAS_microbiome/gemma/gemmaWGS/ids_bothWGS.txt')
nonwgs_ids<-rownames(ps@sam_data)[which(!(rownames(ps@sam_data) %in% wgs_ids$V1))]
points(ps@sam_data$Longitude[which(rownames(ps@sam_data) %in% nonwgs_ids)], 
       ps@sam_data$Latitude[which(rownames(ps@sam_data) %in% nonwgs_ids)], pch=21, bg='black',col='white', lwd=.7, cex=1.7)
points(ps@sam_data$Longitude[which(rownames(ps@sam_data) %in% wgs_ids$V1)], 
       ps@sam_data$Latitude[which(rownames(ps@sam_data) %in% wgs_ids$V1)], pch=23, bg='darkgray',col='black', lwd=.7, cex=1)
map.scale(cex=0.26, ratio=F, x=-122, y=56)
dev.off()

#---- two ASVs for gemma analyses
#make df of asv counts to add to .fam file for gemma
asvs_for_gemma3<-psq@otu_table
asvs_for_gemma3_pheno<-data.frame(id=rownames(psq@sam_data))
identical(colnames(asvs_for_gemma3),asvs_for_gemma3_pheno$id) #TRUE
#beij ASV
asvs_for_gemma3_pheno$beij=as.vector(asvs_for_gemma3['9477d45e83a2f8d1f61fd3cb9db73a0f',])
#micro ASV
asvs_for_gemma3_pheno$micro=as.vector(asvs_for_gemma3['453cf07d3e96e304f23fbde7b355fe70',])

#make sure ids are in correct order
ids_both2<-read.table('../../gemma/ids_both2.txt')
ids_both2$V1<-substr(ids_both2$V1, 1,7)
identical(asvs_for_gemma3_pheno$id,ids_both2$V1) #TRUE

#write files with space sep (beij will be column 16, micro is column 17)
write.table(asvs_for_gemma3_pheno[,c(2:3)], 'asvs_for_gemma3_pheno.txt', sep=' ', quote=F, row.names=F, col.names=F)

#make binary
asvs_for_gemma4_pheno<-asvs_for_gemma3_pheno
asvs_for_gemma4_pheno$beij[which(asvs_for_gemma4_pheno$beij>0)]<-1
asvs_for_gemma4_pheno$micro[which(asvs_for_gemma4_pheno$micro>0)]<-1
write.table(asvs_for_gemma4_pheno[,c(2:3)], '../../gemma/gemma4/asvs_for_gemma4_pheno.txt', sep=' ', quote=F, row.names=F, col.names=F)

#in CCR: 
#paste -d ' ' in_file_gemma.fam asvs_for_gemma3_pheno.txt > test
#mv test in_file_gemma.fam


#----- make plot of alpha ~ year
png('figs/supp/alpha_year.png', width=7, height=7, units='in', res=2400)
par(mfrow=c(2,2), mar=c(3,4,3,1))

boxplot(alpha_sex$Chao1~alpha_sex$Year,xlab='', ylab='Microbiome alpha diversity (Chao1 index)',
        main='', boxwex=0.5, outline=F, ylim=c(0,400), cex.axis=0.8)
stripchart(alpha_sex$Chao1[which(alpha_sex$Year=='2017')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T)
stripchart(alpha_sex$Chao1[which(alpha_sex$Year=='2018')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=2)
stripchart(alpha_sex$Chao1[which(alpha_sex$Year=='2019')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=3)
stripchart(alpha_sex$Chao1[which(alpha_sex$Year=='2020')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=4)
stripchart(alpha_sex$Chao1[which(alpha_sex$Year=='2021')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=5)
stripchart(alpha_sex$Chao1[which(alpha_sex$Year=='2022')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=6)

boxplot(alpha_sex$Shannon~alpha_sex$Year,xlab='', ylab='Microbiome alpha diversity (Shannon)',
        main='', boxwex=0.5, outline=F, ylim=c(0,5.5), cex.axis=0.8)
stripchart(alpha_sex$Shannon[which(alpha_sex$Year=='2017')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T)
stripchart(alpha_sex$Shannon[which(alpha_sex$Year=='2018')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=2)
stripchart(alpha_sex$Shannon[which(alpha_sex$Year=='2019')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=3)
stripchart(alpha_sex$Shannon[which(alpha_sex$Year=='2020')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=4)
stripchart(alpha_sex$Shannon[which(alpha_sex$Year=='2021')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=5)
stripchart(alpha_sex$Shannon[which(alpha_sex$Year=='2022')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=6)

boxplot(alpha_sex$PD~alpha_sex$Year,xlab='', ylab="Microbiome alpha diversity (Faith's PD)",
        main='', boxwex=0.5, outline=F, ylim=c(0,27), cex.axis=0.8)
stripchart(alpha_sex$PD[which(alpha_sex$Year=='2017')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T)
stripchart(alpha_sex$PD[which(alpha_sex$Year=='2018')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=2)
stripchart(alpha_sex$PD[which(alpha_sex$Year=='2019')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=3)
stripchart(alpha_sex$PD[which(alpha_sex$Year=='2020')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=4)
stripchart(alpha_sex$PD[which(alpha_sex$Year=='2021')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=5)
stripchart(alpha_sex$PD[which(alpha_sex$Year=='2022')],vertical=T,pch=21,bg=alpha('gray',0.7), method='jitter',add=T, at=6)

dev.off()



save.image('phyloseq_yrwa.RData')
