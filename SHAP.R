#SHAP project
#Data analysis

#Created by Priscilla Sinclair and Laura Rolon
#Last updated: 01/25/2022  MLR

#Attach ibrarries
library(ggplot2)
library(dplyr)
library(readxl)
library(psych)
library(agricolae)
library(viridis)
library(cowplot)

#Set working directory
setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/Priscilla/quaisbio-main")

#### Spot-on-lawn assay ####
#Import data
spot<-read_xlsx("Inhibition results.xlsx", sheet=1, col_names = TRUE)

#Subset by protective culture
spot_LL<-subset(spot, PC=="LL")
spot_ED<-subset(spot, PC=="ED")

#Add column specifying combinations of temperature and concentration
spot_LL$factorAB <- with(spot_LL, interaction(Concentration, Temp))
spot_ED$factorAB <- with(spot_ED, interaction(Concentration, Temp))

#Calculate statistics
LL_inhib<-describeBy(spot_LL$`Inhibition(mm)`, group=list(spot_LL$Temp, spot_LL$Concentration), mat = TRUE)
ED_inhib<-describeBy(spot_ED$`Inhibition(mm)`, group=list(spot_ED$Temp, spot_ED$Concentration), mat = TRUE)


#ANOVA for the interaction effect (since I care about each independent variable)
anova_LL_inhib<-aov(`Inhibition(mm)` ~ factorAB, data=spot_LL)
summary(anova_LL_inhib)

anova_ED_inhib<-aov(`Inhibition(mm)` ~ factorAB, data=spot_ED)
summary(anova_ED_inhib)

#tukey test
tukey_LL_inhib<-HSD.test(anova_LL_inhib, trt="factorAB") 
tukey_LL_inhib

tukey_ED_inhib<-HSD.test(anova_ED_inhib, trt="factorAB") 
tukey_ED_inhib

#Make dataframe with Tukey groups and order in the way as stat data
tukey_LL_inhib_groups<-tukey_LL_inhib$groups
tukey_LL_inhib_groups$Sample<-as.character(rownames(tukey_LL_inhib_groups))
tukey_LL_inhib_groups<-tukey_LL_inhib_groups[with(tukey_LL_inhib_groups, order(Sample)),]
LL_inhib$TukeyGroups<-tukey_LL_inhib_groups$groups

tukey_ED_inhib_groups<-tukey_ED_inhib$groups
tukey_ED_inhib_groups$Sample<-as.character(rownames(tukey_ED_inhib_groups))
tukey_ED_inhib_groups<-tukey_ED_inhib_groups[with(tukey_ED_inhib_groups, order(Sample)),]
ED_inhib$TukeyGroups<-tukey_ED_inhib_groups$groups

#Plots
spot_LL_plot<- ggplot(LL_inhib, aes(x = group1, y = mean , fill = group2))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("Inhibition (mm)")+xlab("Temperature (??C)")+labs(fill="Lawn concentration")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("LL spot inhibition")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,5) ,breaks= seq(0,5, 0.5))+
  scale_fill_manual(values = c("#BEC8E6", "#8190C8"))
spot_LL_plot
ggsave("spotLL.png", plot=spot_LL_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("spotLL.svg", plot=spot_LL_plot, device="svg", width=5, height=4, units="in", dpi=600)


spot_ED_plot<- ggplot(ED_inhib, aes(x = group1, y = mean , fill = group2))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("Inhibition (mm)")+xlab("Temperature (??C)")+labs(fill="Lawn concentration")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fiED = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill = "transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("ED spot inhibition")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,5) ,breaks= seq(0,5, 0.5))+
  scale_fill_manual(values = c("#FDBF6F", "#E9842B"))
spot_ED_plot
ggsave("spotED.png", plot=spot_ED_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("spotED.svg", plot=spot_ED_plot, device="svg", width=5, height=4, units="in", dpi=600)

#### Supernatant tests ####
#Import data
supernat<-read_xlsx("Inhibition results.xlsx", sheet=2, col_names = TRUE)

#Subset by protective culture
supernat_LL<-subset(supernat, PC=="LL")
supernat_ED<-subset(supernat, PC=="Ed")

#Add column specifying combinations of temperature and concentration
supernat_LL$factorAB <- with(supernat_LL, interaction(Treatment,Time))
supernat_ED$factorAB <- with(supernat_ED, interaction(Treatment,Time))

#Calculate statistics
LL_supernat<-describeBy(supernat_LL$`Inhibition(mm)`, group=list(supernat_LL$Time, supernat_LL$Treatment), mat = TRUE)
ED_supernat<-describeBy(supernat_ED$`Inhibition(mm)`, group=list(supernat_ED$Time, supernat_ED$Treatment), mat = TRUE)


#ANOVA for the interaction effect (since I care about each independent variable)
anova_LL_supernat<-aov(`Inhibition(mm)` ~ factorAB, data=supernat_LL)
summary(anova_LL_supernat)

anova_ED_supernat<-aov(`Inhibition(mm)` ~ factorAB, data=supernat_ED)
summary(anova_ED_supernat)

#tukey test
tukey_LL_supernat<-HSD.test(anova_LL_supernat, trt="factorAB") 
tukey_LL_supernat

tukey_ED_supernat<-HSD.test(anova_ED_supernat, trt="factorAB") 
tukey_ED_supernat

#Make dataframe with Tukey groups and order in the way as stat data
tukey_LL_supernat_groups<-tukey_LL_supernat$groups
tukey_LL_supernat_groups$Sample<-as.character(rownames(tukey_LL_supernat_groups))
tukey_LL_supernat_groups<-tukey_LL_supernat_groups[with(tukey_LL_supernat_groups, order(Sample)),]
LL_supernat$TukeyGroups<-tukey_LL_supernat_groups$groups
LL_supernat$SampleOrder<-c(rep(3,2),rep(1,2), rep(2,2), rep(4,2))
LL_supernat$group2<-as.factor(LL_supernat$group2)
levels(LL_supernat$group2) <- gsub(" ", "\n", levels(LL_supernat$group2))

tukey_ED_supernat_groups<-tukey_ED_supernat$groups
tukey_ED_supernat_groups$Sample<-as.character(rownames(tukey_ED_supernat_groups))
tukey_ED_supernat_groups<-tukey_ED_supernat_groups[with(tukey_ED_supernat_groups, order(Sample)),]
ED_supernat$TukeyGroups<-tukey_ED_supernat_groups$groups
ED_supernat$SampleOrder<-c(rep(3,2),rep(1,2), rep(2,2), rep(4,2))
ED_supernat$group2<-as.factor(ED_supernat$group2)
levels(ED_supernat$group2) <- gsub(" ", "\n", levels(ED_supernat$group2))

#Plots
supernat_LL_plot<- ggplot(LL_supernat, aes(x = reorder(group2,SampleOrder), y = mean , fill = group1))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("Inhibition (mm)")+xlab("Treatment")+labs(fill="Incubation time")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("LL supernat inhibition")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,5) ,breaks= seq(0,5, 0.5))+
  scale_fill_manual(values = c("#8190C8", "#BEC8E6"))
supernat_LL_plot
ggsave("supernatLL.png", plot=supernat_LL_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("supernatLL.svg", plot=supernat_LL_plot, device="svg", width=5, height=4, units="in", dpi=600)


supernat_ED_plot<- ggplot(ED_supernat, aes(x = reorder(group2,SampleOrder), y = mean , fill = group1))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("Inhibition (mm)")+xlab("Treatment")+labs(fill="Incubation time")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("ED supernat inhibition")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_y_continuous(limits=c(0,5) ,breaks= seq(0,5, 0.5))+
  scale_fill_manual(values = c("#E9842B", "#FDBF6F"))
supernat_ED_plot
ggsave("supernatED.png", plot=supernat_ED_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("supernatED.svg", plot=supernat_ED_plot, device="svg", width=5, height=4, units="in", dpi=600)



#Combine plots into fig 1
fig1<-plot_grid(spot_LL_plot, supernat_LL_plot, spot_ED_plot, supernat_ED_plot,
          nrow=2, ncol=2, labels=c("A","B","C","D"), label_size=20)
fig1

ggsave("fig1.png", plot=fig1, device="png", width=9, height=8, units="in", dpi=600)
ggsave("fig1.svg", plot=fig1, device="svg", width=9, height=8, units="in", dpi=600)



#### Biofilm experiments ####

#Import data
biofilm<-read_xlsx("Biofilm results.xlsx", sheet=1, col_names = TRUE, na="NA")

#Subset by endpoint of experiments
biofilm3<-subset(biofilm, Endpoint==3)
biofilm5<-subset(biofilm, Endpoint==5)
biofilm15<-subset(biofilm, Endpoint==15)

#Add column specifying combinations of facility and treatment
biofilm3$factorAB <- with(biofilm3, interaction(Facility, Treatment))
biofilm5$factorAB <- with(biofilm5, interaction(Facility,Treatment))
biofilm15$factorAB <- with(biofilm15, interaction(Facility,Treatment))

#APC data analysis
#Calculate statistics
biofilm3_apc<-describeBy(biofilm3$APC, group=list( biofilm3$Treatment, biofilm3$Facility), mat = TRUE)
biofilm5_apc<-describeBy(biofilm5$APC, group=list(biofilm5$Treatment, biofilm5$Facility), mat = TRUE)
biofilm15_apc<-describeBy(biofilm15$APC, group=list(biofilm15$Treatment, biofilm15$Facility), mat = TRUE)

#ANOVA for the interaction effect (since I care about each independent variable)
anova_biofilm3_apc<-aov(APC ~ factorAB, data=biofilm3)
summary(anova_biofilm3_apc)

anova_biofilm5_apc<-aov(APC ~ factorAB, data=biofilm5)
summary(anova_biofilm5_apc)

anova_biofilm15_apc<-aov(APC ~ factorAB, data=biofilm15)
summary(anova_biofilm15_apc)

#tukey test
tukey_biofilm3_apc<-HSD.test(anova_biofilm3_apc, trt="factorAB") 
tukey_biofilm3_apc

tukey_biofilm5_apc<-HSD.test(anova_biofilm5_apc, trt="factorAB") 
tukey_biofilm5_apc

tukey_biofilm15_apc<-HSD.test(anova_biofilm15_apc, trt="factorAB") 
tukey_biofilm15_apc

#Make dataframe with Tukey groups and order in the way as stat data
tukey_biofilm3_apc_groups<-tukey_biofilm3_apc$groups
tukey_biofilm3_apc_groups$Sample<-as.character(rownames(tukey_biofilm3_apc_groups))
tukey_biofilm3_apc_groups<-tukey_biofilm3_apc_groups[with(tukey_biofilm3_apc_groups, order(Sample)),]
biofilm3_apc$TukeyGroups<-tukey_biofilm3_apc_groups$groups
biofilm3_apc$SampleOrder<-rep(c(4,3,1,2),3)

tukey_biofilm5_apc_groups<-tukey_biofilm5_apc$groups
tukey_biofilm5_apc_groups$Sample<-as.character(rownames(tukey_biofilm5_apc_groups))
tukey_biofilm5_apc_groups<-tukey_biofilm5_apc_groups[with(tukey_biofilm5_apc_groups, order(Sample)),]
biofilm5_apc$TukeyGroups<-tukey_biofilm5_apc_groups$groups
biofilm5_apc$SampleOrder<-rep(c(4,3,1,2),3)

tukey_biofilm15_apc_groups<-tukey_biofilm15_apc$groups
tukey_biofilm15_apc_groups$Sample<-as.character(rownames(tukey_biofilm15_apc_groups))
tukey_biofilm15_apc_groups<-tukey_biofilm15_apc_groups[with(tukey_biofilm15_apc_groups, order(Sample)),]
biofilm15_apc$TukeyGroups<-tukey_biofilm15_apc_groups$groups
biofilm15_apc$SampleOrder<-rep(c(4,3,1,2),3)

#Plots
biofilm3_apc_plot<- ggplot(biofilm3_apc, aes(x = reorder(group1,SampleOrder), y = mean , fill = group2))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + facet_grid(group2~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("log CFU/ml")+xlab("Treatment")+labs(fill="Facility")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Biofilm study - 3 day", subtitle="Aerobic plate count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15), plot.subtitle = element_text(face = 'italic', size=10))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_manual(values = c("#8190C8", "#F3766E", "#E9842B"))
biofilm3_apc_plot
ggsave("biofilm3_apc.png", plot=biofilm3_apc_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("biofilm3_apc.svg", plot=biofilm3_apc_plot, device="svg", width=5, height=4, units="in", dpi=600)

biofilm5_apc_plot<- ggplot(biofilm5_apc, aes(x = reorder(group1,SampleOrder), y = mean , fill = group2))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + facet_grid(group2~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
    ylab("log CFU/ml")+xlab("Treatment")+labs(fill="Facility")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Biofilm study - 5 day", subtitle="Aerobic plate count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15), plot.subtitle = element_text(face = 'italic', size=10))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_manual(values = c("#8190C8", "#F3766E", "#E9842B"))
biofilm5_apc_plot
ggsave("biofilm5_apc.png", plot=biofilm5_apc_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("biofilm5_apc.svg", plot=biofilm5_apc_plot, device="svg", width=5, height=4, units="in", dpi=600)

biofilm15_apc_plot<- ggplot(biofilm15_apc, aes(x = reorder(group1,SampleOrder), y = mean , fill = group2))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + facet_grid(group2~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  ylab("log CFU/ml")+xlab("Treatment")+labs(fill="Facility")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Biofilm study - 15 day", subtitle="Aerobic plate count")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15), plot.subtitle = element_text(face = 'italic', size=10))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_manual(values = c("#8190C8", "#F3766E", "#E9842B"))
biofilm15_apc_plot
ggsave("biofilm15_apc.png", plot=biofilm15_apc_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("biofilm15_apc.svg", plot=biofilm15_apc_plot, device="svg", width=5, height=4, units="in", dpi=600)



#MPN data analysis
#Calculate statistics
biofilm3_lm<-describeBy(biofilm3$LmMPN, group=list( biofilm3$Treatment, biofilm3$Facility), mat = TRUE)
biofilm5_lm<-describeBy(biofilm5$LmMPN, group=list(biofilm5$Treatment, biofilm5$Facility), mat = TRUE)
biofilm15_lm<-describeBy(biofilm15$LmMPN, group=list(biofilm15$Treatment, biofilm15$Facility), mat = TRUE)

#ANOVA for the interaction effect (since I care about each independent variable)
anova_biofilm3_lm<-aov(LmMPN ~ factorAB, data=biofilm3)
summary(anova_biofilm3_lm)

anova_biofilm5_lm<-aov(LmMPN ~ factorAB, data=biofilm5)
summary(anova_biofilm5_lm)

anova_biofilm15_lm<-aov(LmMPN ~ factorAB, data=biofilm15)
summary(anova_biofilm15_lm)

#tukey test
tukey_biofilm3_lm<-HSD.test(anova_biofilm3_lm, trt="factorAB") 
tukey_biofilm3_lm

tukey_biofilm5_lm<-HSD.test(anova_biofilm5_lm, trt="factorAB") 
tukey_biofilm5_lm

tukey_biofilm15_lm<-HSD.test(anova_biofilm15_lm, trt="factorAB") 
tukey_biofilm15_lm

#Make dataframe with Tukey groups and order in the way as stat data
tukey_biofilm3_lm_groups<-tukey_biofilm3_lm$groups
tukey_biofilm3_lm_groups$Sample<-as.character(rownames(tukey_biofilm3_lm_groups))
tukey_biofilm3_lm_groups<-tukey_biofilm3_lm_groups[with(tukey_biofilm3_lm_groups, order(Sample)),]
biofilm3_lm$TukeyGroups<-tukey_biofilm3_lm_groups$groups
biofilm3_lm$SampleOrder<-rep(c(4,3,1,2),3)

tukey_biofilm5_lm_groups<-tukey_biofilm5_lm$groups
tukey_biofilm5_lm_groups$Sample<-as.character(rownames(tukey_biofilm5_lm_groups))
tukey_biofilm5_lm_groups<-tukey_biofilm5_lm_groups[with(tukey_biofilm5_lm_groups, order(Sample)),]
biofilm5_lm$TukeyGroups<-tukey_biofilm5_lm_groups$groups
biofilm5_lm$SampleOrder<-rep(c(4,3,1,2),3)

tukey_biofilm15_lm_groups<-tukey_biofilm15_lm$groups
tukey_biofilm15_lm_groups$Sample<-as.character(rownames(tukey_biofilm15_lm_groups))
tukey_biofilm15_lm_groups<-tukey_biofilm15_lm_groups[with(tukey_biofilm15_lm_groups, order(Sample)),]
biofilm15_lm$TukeyGroups<-tukey_biofilm15_lm_groups$groups
biofilm15_lm$SampleOrder<-rep(c(4,3,1,2),3)

#Plots
biofilm3_lm_plot<- ggplot(biofilm3_lm, aes(x = reorder(group1,SampleOrder), y = mean , fill = group2))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + facet_grid(group2~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  geom_hline(yintercept=1.522, color='grey20', linetype='dashed')+
  ylab("log MPN/ml")+xlab("Treatment")+labs(fill="Facility")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Biofilm study - 3 day", subtitle="L. monocytogenes quantification")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15), plot.subtitle = element_text(face = 'italic', size=10))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_manual(values = c("#8190C8", "#F3766E", "#E9842B"))
biofilm3_lm_plot
ggsave("biofilm3_lm.png", plot=biofilm3_lm_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("biofilm3_lm.svg", plot=biofilm3_lm_plot, device="svg", width=5, height=4, units="in", dpi=600)

biofilm5_lm_plot<- ggplot(biofilm5_lm, aes(x = reorder(group1,SampleOrder), y = mean , fill = group2))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + facet_grid(group2~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  geom_hline(yintercept=1.522, color='grey20', linetype='dashed')+
  ylab("log MPN/ml")+xlab("Treatment")+labs(fill="Facility")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Biofilm study - 5 day", subtitle="L. monocytogenes quantification")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15), plot.subtitle = element_text(face = 'italic', size=10))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_manual(values = c("#8190C8", "#F3766E", "#E9842B"))
biofilm5_lm_plot
ggsave("biofilm5_lm.png", plot=biofilm5_lm_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("biofilm5_lm.svg", plot=biofilm5_lm_plot, device="svg", width=5, height=4, units="in", dpi=600)

biofilm15_lm_plot<- ggplot(biofilm15_lm, aes(x = reorder(group1,SampleOrder), y = mean , fill = group2))  + 
  geom_bar(stat = "identity", color='black', position = position_dodge()) + facet_grid(group2~.)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.1,position=position_dodge(.9))+
  geom_text(aes(label=TukeyGroups), position=position_dodge(width=0.9), vjust=-0.25)+
  geom_hline(yintercept=1.522, color='grey20', linetype='dashed')+
  ylab("log MPN/ml")+xlab("Treatment")+labs(fill="Facility")+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9), legend.position = 'right') +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 1)) +
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Biofilm study - 15 day", subtitle="L. monocytogenes quantification")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15), plot.subtitle = element_text(face = 'italic', size=10))+
  scale_y_continuous(limits=c(0,10) ,breaks= seq(0,10, 2))+
  scale_fill_manual(values = c("#8190C8", "#F3766E", "#E9842B"))
biofilm15_lm_plot
ggsave("biofilm15_lm.png", plot=biofilm15_lm_plot, device="png", width=5, height=4, units="in", dpi=600)
ggsave("biofilm15_lm.svg", plot=biofilm15_lm_plot, device="svg", width=5, height=4, units="in", dpi=600)


#Combine figures
fig2<-plot_grid(biofilm3_apc_plot, biofilm3_lm_plot, biofilm5_apc_plot, biofilm5_lm_plot, biofilm15_apc_plot ,biofilm15_lm_plot,
                nrow=3, ncol=2, labels=c("A","B","C","D","E","F"), label_size=20)
fig2

ggsave("fig2.png", plot=fig2, device="png", width=10, height=12, units="in", dpi=600)
ggsave("fig2.svg", plot=fig2, device="svg", width=10, height=12, units="in", dpi=600)