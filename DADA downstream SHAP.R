#SHAP microbiome data analysis 
#Got ASVs from DADA2 pipeline using default parameters.

#Last updated: MLR 1/25/2022

#Set working directory
setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/Priscilla/quaisbio-main/Microbiome/")

#Attach libraries
library(ggplot2)
library(dplyr)
library(phyloseq)
library(zCompositions)
library(compositions)
library(viridis)
library(svglite)
library(pairwiseAdonis)
library(decontam)
library(readxl)

#### Import data ####
asvs_all<-read.csv('ASV.csv', header = TRUE, row.names = 1)
taxon_all<-read.csv('Taxon.csv', header = TRUE, row.names = 1)
metadata_all<-as.data.frame(read_excel('Biofilm metadata.xlsx', col_names = TRUE, sheet=1))
rownames(metadata_all)<-metadata_all$SampleID

#Add '_unclassified' marker to NAs in taxon table
taxon_all$Phylum<-ifelse(is.na(taxon_all$Phylum), paste(taxon_all$Kingdom, "unclassified", sep = '_'), taxon_all$Phylum)
taxon_all$Class<-ifelse(is.na(taxon_all$Class), paste(taxon_all$Phylum, "unclassified", sep = '_'), taxon_all$Class)
taxon_all$Order<-ifelse(is.na(taxon_all$Order), paste(taxon_all$Class, "unclassified", sep = '_'), taxon_all$Order)
taxon_all$Family<-ifelse(is.na(taxon_all$Family), paste(taxon_all$Order, "unclassified", sep = '_'), taxon_all$Family)
taxon_all$Genus<-ifelse(is.na(taxon_all$Genus), paste(taxon_all$Family, "unclassified", sep = '_'), taxon_all$Genus)
taxon_all$Species<-ifelse(is.na(taxon_all$Species), paste(taxon_all$Genus, "unclassified", sep = '_'), taxon_all$Species)

#Remove extra _unclassified
taxon_all$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Class)
taxon_all$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Order)
taxon_all$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Order)
taxon_all$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Species)

#Convert asv and taxon tables to matrix
asvs_all<-as.matrix(asvs_all)
taxon_all<-as.matrix(taxon_all)

#Make phylose object
ps_all<-phyloseq(otu_table(asvs_all, taxa_are_rows = FALSE), tax_table(taxon_all), sample_data(metadata_all))
# 
# #Split phyloseq by 
# ps_Obj1<-subset_samples(ps_all, Obj == 1)
# 
# #Get ASV ad metadata table from phyloseq object
# asv_obj1<-as.data.frame(t(otu_table(ps_Obj1)))
# meta_obj1<-sample_data(ps_Obj1)
# 
# #Remove ASVs with zero counts in all samples
# asv_obj1<-asv_obj1[ which(rowSums(asv_obj1)>0),]


# # RAREFACTION CURVES
# #Add ggrare function to R. Code available at https://rdrr.io/github/gauravsk/ranacapa/src/R/ggrare.R
# phyloseq_rare_obj1 = phyloseq(otu_table(asv_obj1, taxa_are_rows = TRUE), tax_table(taxon_all), sample_data(meta_obj1))
# rare_obj1 <- ggrare(phyloseq_rare_obj1, step = 100, se=TRUE, color="Facility")
# 
# #Plot rarefaction curves
# rare_plot_obj1 <- rare_obj1 + 
#   theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
#   #xlab("Number of reads (x 1,000)") + ylab("Number of unique OTUs") +
#   scale_x_continuous(breaks= seq(0,100000, 5000), labels = function(x){x/1000}) + ylim(0,500)+
#   theme(axis.text = element_text(color='black', size=8)) +
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA))+
#   theme(panel.grid.major = element_line(color='gray'))+
#   theme(legend.position = "bottom")+
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
#   scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
#   scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
# rare_plot_obj1
# ggsave("Rarefaction ASV Obj1.png", plot =rare_plot_obj1, device="png", width=5, height=4, units="in",dpi=600)
# ggsave("Rarefaction ASV Obj1.svg", plot =rare_plot_obj1, device="svg", width=5, height=4, units="in",dpi=600)
# 
# 
# # Estimate % discovered diversity 
# #Calculate estimated richness for each sample.
# ## Original script adopted from "https://cran.r-project.org/web/packages/SpadeR/SpadeR.pdf"
# ## SpadeR::ChaoSpecies function is used to estimate richness.
# 
# richness_estimate = function(otu,...) {
#   options(warn = -1)
#   b = data.frame(matrix(nrow=as.matrix(dim(otu))[2,], ncol=3))
#   colnames(b) <- c("Chao1 Estimates","Observed OTUs", "%Covered Species")
#   for (i in 1:as.matrix(dim(otu))[2,]) {
#     a =SpadeR::ChaoSpecies(otu[,i], datatype="abundance", k=10, conf=0.95)
#     b[i,1]= as.numeric(a$Species_table[3,1])
#     b[i,2]= apply(as.data.frame(a$Basic_data_information),2,as.numeric)[2,2]
#     b[i,3]= (b[i,2]/b[i,1])*100
#     rownames(b) <- colnames(otu) }
#   print(b)
# }
# 
# #Estimate richness using Chao1 index and % discovery
# 
# #Add Taejung if I need to adapt the code for working with an ASV table
# 
# asv_obj1<-as.data.frame(asv_obj1)
# 
# SpadeR::ChaoSpecies(asv_obj1, datatype="abundance", k=10, conf=0.95)
# 
# spadeR_estimate_obj1 <- richness_estimate(asv_obj1)
# spadeR_estimate_obj1
# 
# write.csv(spadeR_estimate_obj1, file='discovered_diversity Obj1.csv')

#CHECH CONTROLS
#Note: For this batch of sequences, the NC did not amplify so it was not included in the sequencing experiment.

#Subset Controls
phyloseq_PC<-subset_samples(ps_all, SampleID=="PC")

#Make long table
phyloseq_PC<-psmelt(phyloseq_PC)

#Remove rows with less than 100 reads
phyloseq_PC<-subset(phyloseq_PC, Abundance>100)

#Plot reads of control by PC, PC_DNA, NC
barplot_PC<-ggplot(phyloseq_PC, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=Facility))+
  geom_bar(stat='identity', color='black')+ theme()+
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,65000, 5000), labels = function(x){x/10000}, limits=c(0,65000)) + 
  ggtitle("Obj2 - PC")+ xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_PC.png", plot =barplot_PC, device="png", width=20, height=10, units="in",dpi=600)



#Calculate relative abundance for PC sample
phyloseq_PC$RA<-phyloseq_PC$Abundance*100/sum(phyloseq_PC$Abundance)

phyloseq_PC_over1<-filter(phyloseq_PC, RA>1)

PC_RA<-ggplot(phyloseq_PC_over1, aes(x=SampleID, y=RA, fill=Genus))+
  geom_bar(stat='identity', color='black')
ggsave("Controls_PC_RA.png", plot =PC_RA, device="png", width=5, height=8, units="in",dpi=600)


#Remove positive and negative control from my OTU table
phyloseq_clean<-subset_samples(ps_all, SampleID!="PC")

#Get ASV table from phyloseq object
asv_clean<-as.data.frame(t(otu_table(phyloseq_clean)))
tail(rowSums(asv_clean))

#Remove OTUs with zero counts in all samples
asv_clean<-asv_clean[ which(rowSums(asv_clean)>0),]
asv_clean<-t(asv_clean)

#Get metadata 
metadata_clean<-subset(metadata_all, SampleID!="PC")

#Compositional analysis of microbiome data -based on Microbiome Analysis in R. Chap 10.
#Step 1: Convert OTU table to appropriate format. Following step requires samples on rows and OTUs in columns
head(asv_clean)

#Step 2: Replace zero values before clr transformation. Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0<-t(cmultRepl(asv_clean, label=0, method="CZM", output="p-counts")) #214173 corrected values

head(asv.n0) #output table needs to have samples in columns and OTUs in rows

#Step 3: Convert data to proportions
asv.n0_prop<-apply(asv.n0, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_prop_f<-asv.n0[apply(asv.n0_prop, 1, min) > 0.0000001, ]
head(asv.n0_prop_f) #Check that samples are on columns and asvs in rows

#Step 5: perform CLR transformation
asv.n0.clr<-t(apply(asv.n0, 2, function(x){log(x)-mean(log(x))}))
head(asv.n0.clr) #Check output table. Samples should be in rows and asvs in columns

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr<-prcomp(asv.n0.clr)

png("Screeplot - PCA .png")
par(mar=c(2,2,2,2))
screeplot(pc.clr, type='lines', main="Screeplot - PCA")
dev.off()

#Calculate total variance of the data
mvar.clr<-mvar(asv.n0.clr)

#Display results
row<-rownames(asv.n0.clr) #Make vector with sample names
pc_out<-as.data.frame(pc.clr$x[,1:2]) #Get PC1 and PC2
pc_out_meta<-as.data.frame(bind_cols(pc_out,metadata_clean)) #Add metadata information
row.names(pc_out_meta)<-row #Add rownames to dataframe

# Make PCA plot
PCA_fac <- ggplot(pc_out_meta, aes(x=PC1,y=PC2, color=Facility, shape=Treatment))+
  geom_point(size=5)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13)) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr$sdev[1]^2/mvar.clr*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr$sdev[2]^2/mvar.clr*100, digits=1), "%", sep="")) +
  ggtitle("PCA - by facility and treatment")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_manual(values = c("#8190C8", "#F3766E", "#E9842B"))
PCA_fac
ggsave("PCA_fac.png", plot=PCA_fac, device="png", width=8, height=8, units="in", dpi=600)
ggsave("PCA_fac.svg", plot=PCA_fac, device="svg", width=8, height=8, units="in", dpi=600)

PCA_time<- ggplot(pc_out_meta, aes(x=PC1,y=PC2, color=Endpoint, shape=Treatment))+
  geom_point(size=5)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13)) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr$sdev[1]^2/mvar.clr*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr$sdev[2]^2/mvar.clr*100, digits=1), "%", sep="")) +
  ggtitle("PCA - by facility and treatment")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_manual(values = c("#8190C8", "#F3766E", "#E9842B"))
PCA_time
ggsave("PCA_time.png", plot=PCA_time, device="png", width=8, height=8, units="in", dpi=600)
ggsave("PCA_time.svg", plot=PCA_time, device="svg", width=8, height=8, units="in", dpi=600)

# PERMANOVA
#Calculate Aitchinson distance
dist<-dist(asv.n0.clr, method='euclidean')

#Permanova by endpoint
permanova_end<-pairwise.adonis2(dist~Endpoint, data=metadata_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_end


#Permanova by facility
permanova_fac<-pairwise.adonis2(dist~Facility, data=metadata_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_fac

#Permanova by facility
permanova_trt<-pairwise.adonis2(dist~Treatment, data=metadata_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_trt


##Results: Enpoint is significant, not facility or treatment.

# STACKED BARPLOT 
#Note: used compositional approach to transform the sample counts to compositions. 

#Transform sample counts into compositions
asv.n0.acomp<-as.data.frame(acomp(t(asv.n0)), total=1)
rowSums(asv.n0.acomp) #Verify rows sum to 1

#Make Phyloseq object
phyloseq_RA <- phyloseq(otu_table(asv.n0.acomp, taxa_are_rows = FALSE), tax_table(taxon_all), sample_data(metadata_clean))

#Make long format table from Phyloseq object
asv_long <- phyloseq_RA %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

# Calculate mean relative abundance by Treatment, Facility and Endpoint for each OTU
asv_mean<-asv_long %>%
  group_by(OTU, Treatment, Facility, Endpoint, Family, Genus, SampleOrder)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

#Filter table to obtain only OTUs with over 2% in at least one sample
asv_over2abund <- filter(asv_mean, Mean>2)

#Stacked barplot by treatment at the OTU level
barplot_asv<-ggplot(asv_over2abund, aes(x=reorder(Treatment,SampleOrder), y=Mean, fill=Genus))+
  geom_bar(stat='identity', color='black')+facet_grid(Endpoint~Facility, scales = "free", space = 'free')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=3)+ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 4)) +
  theme(legend.position = "bottom")+ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Microbiota of biofilms")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_microbiome.png", plot=barplot_asv, device="png", width=15, height=13, units="in", dpi=600)
ggsave("Barplots_microbiome.svg", plot=barplot_asv, device="svg", width=15, height=13, units="in", dpi=600)






