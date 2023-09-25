---
title: "Rumen_ARM_ver01"
author: "Gerardo Diaz"
date: "`r Sys.Date()`"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 0. Loading packages

```{r libraries, message=FALSE}
library("phyloseq")
library("forcats")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("microbiome")   # Extra data manipulation
library(metagenomeSeq)  # Data normalization (?) + Statistical testing
library(limma)          # Statistical testing
library("statmod")      # Will work with limma for "duplicateCorrelation" function
library(psych)          # Data exploration
library("tidyverse")    # ggplot2, dplyr, tidyr, forcats, tibble, string, purrr, readr
library(vegan)
library(lme4)           # Statistical testing
library(lmerTest)       # Statistical testing
library(emmeans)        # Statistical testing
library(car)            # Statistical testing
library(table1)         # for doing tables
library(lemon)          # for volcano plot
library("ggrepel")      # for volcano plot

```

# 1. Loading RAW data 

Change your working directory to where the files are located

Three tables are needed

* OTU (ARG)
* Taxonomy (ARG classification)
* Samples

They are read from a single Excel file where each sheet contains one of the tables

WE ARE GOING TO USE THE AMR_analytic_matrix

```{r}
  otu_mat<- read_excel("RumenAMR_phyloseq_Ver02.xlsx", sheet = "OTU_matrix")
  tax_mat<- read_excel("RumenAMR_phyloseq_Ver02.xlsx", sheet = "Taxonomy_table")
  samples_df<- read_excel("RumenAMR_phyloseq_Ver02.xlsx", sheet = "Samples")
    
  samples_df$collection_date= factor(samples_df$collection_date, levels=c("9.21.21", "10.18.21", "10.20.21"))
  samples_df$collection_day= factor(samples_df$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))
   samples_df$castration_group= factor(samples_df$castration_group, levels=c("Birth", "Turnout", "Pre_weaning","Weaning"))
   samples_df$castration= factor(samples_df$castration, levels=c("Birth", "Turnout", "Pre_weaning","Weaning", "Not_castrated"))
   samples_df$castration_status= factor(samples_df$castration_status, levels = c("Castrated", "Not_castrated"))
   samples_df$weaning= factor(samples_df$weaning, levels=c("Fence_line", "Truck","Not_weaned"))
   samples_df$weaning_group=factor(samples_df$weaning_group, levels=c("Fence_line", "Truck"))
   samples_df$weaning_status=factor(samples_df$weaning_status, levels=c("Weaned", "Not_weaned"))
  samples_df$cow_ID= as.factor(samples_df$cow_ID)
  samples_df$extraction_date= as.factor(samples_df$extraction_date)
  samples_df$extraction_run= as.factor(samples_df$extraction_run)

  
```

## 1.1. Formatting data as needed to phyloseq

Phyloseq objects need to have row.names

* define the row names from the otu column in otu_mat and tax_mat

```{r}
  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("ARG") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("ARG")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("samples")
```

Transform into matrixes otu and tax tables (sample table can be left as data frame)

```{r}
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
```

Transform to phyloseq objects

```{r}
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)
  
  Rumen_AMR <- phyloseq(OTU, TAX, samples)
  Rumen_AMR
```  

At the end of this step we have a phyloseq object called **Rumen_AMR**, which contains RAW DATA: 194 ARG classified by 4 classification ranks (Type, Class, Mechanism, Gene group), 99 samples and 32 variables. Since no ARG were found in blank2 (sample R0120), it was removed from analysis.

We will modify (subset, delete ARG, analyze a given Classification, etc) this phyloseq object according to our objectives.

# 2. Quality control of reads

```{r}

#plotting classified & unclassified reads (MEGAres)

reads <- read_excel("../METADATA_ver05.xlsx", sheet = "Reads")
reads <-data.frame(reads)
reads$MEGAres_classification<-factor(reads$MEGAres_classification, levels = c("No classified","Classified"))

ggplot(reads, aes(fill=MEGAres_classification, y=MEGAres_classified, x=Sample)) + 
  geom_bar(position="fill", stat="identity")+theme_bw()+facet_wrap(~Collection_day, nrow=1, scales="free_x")+scale_fill_manual(values = c("grey", "black"))+
  scale_color_manual(values = c("grey", "black"))+theme_bw()+
  theme(axis.text.x=element_blank()) + labs(title="MEGAres classification", x="Samples", y = "Percentage of reads")
  

# plotting number of reads by weaning, castration and collection day

# taking only samples and "as.factor"
samples_2<-subset(samples_df, samples_df$sample_type=="Sample")
samples_2=data.frame(samples_2)
samples_2$castration_group=factor(samples_2$castration_group, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning"))
samples_2$collection_day=factor(samples_2$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))

samples_D1<-subset(samples_2, samples_2$collection_day=="Pre_weaning")
samples_D1=data.frame(samples_D1)
samples_D2<-subset(samples_2, samples_2$collection_day=="At_weaning")
samples_D2=data.frame(samples_D2)
samples_D3<-subset(samples_2, samples_2$collection_day=="Post_weaning")
samples_D3=data.frame(samples_D3)


# Plotting
   ggplot(samples_2, aes(castration_group,seq_AMRclassifiedRAW,        
                         fill=castration_group))+
                 geom_boxplot() +facet_wrap(~collection_day, nrow=1)+
   scale_fill_manual(values = c("#FF3333", "blue", "#F3F304", "#009900"))+
  scale_color_manual(values = c("#FF3333", "blue", "#F3F304", "#009900"))+
   theme_bw()+ geom_jitter(width=0.001)+
   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
   theme(axis.text.x=element_blank())+
   ylab("Classified reads (hits)")+xlab("Collection day")



```


# 3. Cleaning/managing dataset according to your preference

## 3.1. Subsetting phyloseq object for SAMPLES only

Subsetting controls:
- POS-CONTORL: Sub-setting the mock community. Phyloseq object obtained: **Rumen_Mock**
- NEG-CONTROL: Sub-setting the DNA extraction blanks. Phyloseq object obtained: **Rumen_Blank**
- SAMPLES: Sub-setting only samples (95). Phyloseq object obtained: **Rumen_AMR_Sample**

```{r}
  Rumen_Mock <- subset_samples(Rumen_AMR, sample_type =="Mock")
  Rumen_Blank <- subset_samples(Rumen_AMR,sample_type =="Blank") 
  Rumen_AMR_Sample <- subset_samples(Rumen_AMR, sample_type =="Sample")
  
```

**Rumen_Mock** will have 2 samples
**Rumen_Blank** will have 2 samples (Sample R0120 was deleted because no ARG)
**Rumen_AMR_Sample** will have 95 samples

We can also subset by SAMPLING DATE:
- 1st day: **Rumen_1**
- 2nd day: **Rumen_2**
- 3rd day: **Rumen_3**

```{r}
  
  Rumen_1 <- subset_samples(Rumen_AMR_Sample, collection_day =="Pre_weaning")
  Rumen_2 <- subset_samples(Rumen_AMR_Sample, collection_day =="At_weaning")
  Rumen_3 <- subset_samples(Rumen_AMR_Sample, collection_day =="Post_weaning")
 
```

**Rumen_1** will have 32 samples
**Rumen_2** will have 32 samples
**Rumen_3** will have 31 samples

## 3.2. Normalizing data

### 3.2.1. Cumulative sum scaling (CSS)

You can do it using metagenomeSeq package
We will use the phyloseq object **Rumen_AMR_Sample**

```{r}

# Changing from phyloseq object to metagenomeSeq object for CSS norm
rumen.metaseq <- phyloseq_to_metagenomeSeq(Rumen_AMR_Sample) # change here

# CSS normalization
rumen.metaseq.norm<- cumNorm(rumen.metaseq, p=cumNormStat(rumen.metaseq))
CSS_rumen.metaseq <- MRcounts(rumen.metaseq.norm, norm = TRUE)

# Changing back from metagenomeSeq object to phyloseq object
CSS_rumen <- merge_phyloseq(otu_table(CSS_rumen.metaseq,
taxa_are_rows=T),sample_data(Rumen_AMR_Sample),tax_table(Rumen_AMR_Sample))

# show your CSS normalized phyloseq object
head(otu_table(CSS_rumen))
CSS_rumen

```


## 3.3. Agglomerate ARG to different ranks

Agglomerate taxa to **Gene group level**

```{r}
  Rumen_Gene= tax_glom(CSS_rumen, "Gene_group")
  Rumen_Gene
```

You may get a **Rumen_Gene** with 111 taxa and 95 samples

# 4. RELATIVE ABUNDANCE (ARG barplots)

## 4.1.  Plotting the Gene group ARG divided by weaning

```{r}


# 1. Converting to Percentage and transforming Genus tax_glom object to long-format table
  phy_relative3 <- transform_sample_counts(Rumen_Gene, function(x) (x / sum(x))*100 )
  phy_relative_long3 <- psmelt(phy_relative3)

# 2. Creating a new column (variable) for the  long-format table, called "mean_relative_abund" 
  phy_relative_long3 <- phy_relative_long3 %>%
  group_by(Gene_group) %>%
  mutate(mean_relative_abund = mean(Abundance))

# 3. Formatting the "Genus" column variable as names and the "mean_relative_abund" column variable as continuous numbers
  phy_relative_long3$Gene_group <- as.character(phy_relative_long3$Gene_group)
  phy_relative_long3$mean_relative_abund <- as.numeric(phy_relative_long3$mean_relative_abund)

# 4. Giving the name "Others (< 1%)" to all the OTUs in the "Genus" column that have less than 1 in the "Abundace" column in the long-format table 
phy_relative_long3$Gene_group[phy_relative_long3$Abundance < 1] <- "Others (< 1%)" 

# 5. PLOTTING by collection day, castration and weaning group


phy_relative_long3 %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Gene_group)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")


```


# 5. ALPHA-DIVERSITY

## 5.1. Obtaining alpha diversity indices (Obs, shannon, Simpson, Evenness)

```{r}
# WITHOUT CSS Normalization, NOT GLOM, Obtain alpha diversity for every sample at the ARG level
genus_alpha.div= estimate_richness(Rumen_AMR, measures=c("Observed", "InvSimpson", "Shannon"))
write.csv(genus_alpha.div, "Alpha_Rumen_AMR.csv")

genus_alpha.div2= evenness(Rumen_AMR, 'pielou')
write.csv(genus_alpha.div2, "Alpha-Evenness_Rumen_AMR.csv")

```

## 5.2. Plotting alpha diversity

```{r}
# Formatting data
samples_2$castration_group=factor(samples_2$castration_group, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning"))
samples_2$castration=factor(samples_2$castration, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning", "Not_castrated"))
samples_2$weaning_group=factor(samples_2$weaning_group, levels=c("Fence_line", "Truck"))
samples_2$weaning=factor(samples_2$weaning, levels=c("Fence_line", "Truck", "Not_weaned"))
samples_2$collection_day=factor(samples_2$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))



# plotting richness, shannon and Pielou by WEANING and Collection day
  ggplot(samples_2, aes(collection_day,alphaAMR_Shan, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ geom_jitter(position=position_dodge(0.8), size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Shannon Index")+xlab("Collection day")
 
  md= samples_2
  md1= samples_D1
md2= samples_D2
md3= samples_D3
  
  
  ModelX = lmer(alphaAMR_Shan ~ weaning+collection_day + castration + age_days + seq_K0.1classified + (1|cow_ID), data=md, REML = T)
  ModelY = lmer(alphaAMR_Shan ~ collection_day+weaning + castration + age_days + seq_K0.1classified + (1|cow_ID), data=md, REML = T)
  ModelX1 = lm(alphaAMR_Shan ~ weaning+castration + age_days + seq_K0.1classified, data=md3)
  
  summary(ModelX)

  anova(ModelX)
  Anova(ModelX, type = "III")
  
  emmeans(ModelX,pairwise~weaning)


  
```

# 6. BETA-DIVERSITY

## 6.1. General

```{r}
# Get a NMDS plot in a bray-curtis matrix for OVERVIEW (GENUS LEVEL)

genus.ord <- ordinate(CSS_rumen, "NMDS", "bray")

plot_ordination(CSS_rumen, genus.ord, type="sample", color="collection_day")+
  geom_point(size=4.5)+scale_fill_manual(values = c("Brown", "Grey", "Black"))+
  scale_color_manual(values = c("Brown", "Grey", "Black"))+theme_classic()

p<-plot_ordination(CSS_rumen, genus.ord, type="sample", color="weaning", shape="collection_day")+
  geom_point(size=4.5, alpha=0.8)+scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+scale_shape_manual(values = c(15,16,17))+  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+theme_classic()

p$layers
p$layers <- p$layers[-1]
p


pp<-plot_ordination(CSS_rumen, genus.ord, type="sample", color="weaning")+
  geom_point(size=4.5, alpha=0.8)+scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+stat_ellipse() +theme_classic()

pp$layers
pp$layers <- pp$layers[-1]
pp


metadata= data.frame(sample_data(CSS_rumen))
dist = vegdist(t(otu_table(CSS_rumen)), method ="bray" )


  set.seed(143)
  permanova5.1 = adonis2(dist ~ weaning + collection_day  ,strata=metadata$collection_day,  data=metadata)  
  permanova5.3 = adonis2(dist ~ collection_day + weaning + castration,strata=metadata$weaning,  data=metadata)


```


# 7. Differential abundance

## 7.1. OBTAINING DATA

```{r}

# 1. DEFINE YOUR METAGENOMESEQ EXPERIMENT (cumNorm-ed)

# we need to filter the OTUs at Genus level to be present in at least in 2 samples (see ?filterData for help). We will use previous MRexperiment used to filter at the Genus level:
  Rumen.metaseq <- phyloseq_to_metagenomeSeq(CSS_rumen) # it will change to MRexperiment as required by the "filterData" function

  Rumen.metaseq

  Rumen_Gen_filtered <- filterData(aggTax(Rumen.metaseq, lvl="Gene_group"), present = 2) # it will aggregate taxa at Genus (and give Genus names instead of only "OTU XXXX"), and filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)

  Rumen_Gen_filtered.metaseq <- cumNorm(Rumen_Gen_filtered) # it will normalize again the new filtered MRexperiment 


# 2. DEFINE YOUR METADATA VARIABLES TO INCLUDE IN THE MODEL (pData)

  collection_day <- pData(Rumen_Gen_filtered.metaseq)$collection_day
  castration <-  pData(Rumen_Gen_filtered.metaseq)$castration # will include 4 castration strategies + Not_castrated
  weaning <-  pData(Rumen_Gen_filtered.metaseq)$weaning # will include 2 weaning strategies + Not_weaned
  cow_ID <- pData(Rumen_Gen_filtered.metaseq)$cow_ID
  
  
# 3. Define your settings (zigControl), we will use the ones by default

  settings = zigControl(maxit = 10, verbose = TRUE)

# 4. DEFINE YOUR MODEL ZERO TO ACCOUNT FOR DEPTH OF COVERAGE (zero_mod)

## MODEL 0: ZERO MODEL FOR ALL DAYS (depth of coverage of the NOT-normalized OTU table from glom phyloseq obj)

  Rumen_RAW_Genus_filtered <- filterData(aggTax(rumen.metaseq,lvl="Gene_group"), present = 2) # it will filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)

  zero_mod_Gen <- model.matrix(~0+log(libSize(Rumen_RAW_Genus_filtered))) # model it!


# 5. DEFINE YOUR MODEL TO TEST (model.matrix), FIT THE MODEL (fitZig), MAKE PAIRWISE CONTRASTS (makeContrasts, contrasts.fit, eBayes) & EXTRACT RESULTS (topTable) 


# MODEL 1: COllection day

# Assign model (model.matrix)
  design_collectionday = model.matrix(~0 + collection_day)
  
# Assign correlated data (random effect)
  v_collectionday_gen <- voom(MRcounts(Rumen_Gen_filtered.metaseq), design_collectionday)

  dup_collectionday_gen <- duplicateCorrelation(v_collectionday_gen, block=cow_ID, design= design_collectionday)

  dup_collectionday_gen$consensus.correlation # to check your correlation coeff?

# Make ZIG model 
  ZIG_Genus_collectionday = fitZig(obj= Rumen_Gen_filtered.metaseq, 
                                   mod = design_collectionday, 
                                   control = settings, 
                                   zeroMod=zero_mod_Gen, useCSSoffset=FALSE,
                                   useMixedModel=dup_collectionday_gen$consensus.correlation, 
                                   block=cow_ID)
  ZIGFit_Genus_collectionday = ZIG_Genus_collectionday@fit
  Mod_Genus_collectionday = ZIG_Genus_collectionday@fit$design

# Contrasts by castration group
  contrast_Genus_collectionday= makeContrasts(collection_dayAt_weaning-collection_dayPre_weaning ,
                                           collection_dayPost_weaning-collection_dayPre_weaning,
                                           collection_dayPost_weaning-collection_dayAt_weaning,
                                           levels=Mod_Genus_collectionday)
  res_Genus_collectionday = contrasts.fit(ZIGFit_Genus_collectionday, contrast_Genus_collectionday)
  resEB_Genus_collectionday = eBayes(res_Genus_collectionday)

# Extract results
  fz_Genus_collectionday <- topTable(resEB_Genus_collectionday, adjust.method="BH",number = 1200)

  fz_Genus_collectionday_At_Pre <- topTable(resEB_Genus_collectionday, adjust.method="BH",coef = 1,number = 1200)
  fz_Genus_collectionday_Post_Pre <- topTable(resEB_Genus_collectionday, adjust.method="BH",coef = 2,number = 1200)
  fz_Genus_collectionday_Post_At <- topTable(resEB_Genus_collectionday, adjust.method="BH",coef = 3,number = 1200)

# Write CSV files
  write.csv(fz_Genus_collectionday_At_Pre, "./fz_Genus_collectionday_At_Pre.csv")
  write.csv(fz_Genus_collectionday_Post_Pre, "./fz_Genus_collectionday_Post_Pre.csv")
  write.csv(fz_Genus_collectionday_Post_At, "./fz_Genus_collectionday_Post_At.csv")
  


# MODEL 3: Weaning overall

# Assign model (model.matrix): same as Phylum
  design_weaning = model.matrix(~0 + weaning)

# Assign correlated data (random effect)
  v_weaning_gen <- voom(MRcounts(Rumen_Gen_filtered.metaseq), design_weaning)

  dup_weaning_gen <- duplicateCorrelation(v_weaning_gen, block=cow_ID, design= design_weaning)

  dup_weaning_gen$consensus.correlation # to check your correlation coeff?

# Make ZIG model 
  ZIG_Genus_weaning = fitZig(obj= Rumen_Gen_filtered.metaseq, 
                             mod = design_weaning, 
                             control = settings, 
                             zeroMod=zero_mod_Gen, useCSSoffset=FALSE, 
                             useMixedModel=dup_weaning_gen$consensus.correlation, 
                             block=cow_ID)
  ZIGFit_Genus_weaning = ZIG_Genus_weaning@fit
  Mod_Genus_weaning = ZIG_Genus_weaning@fit$design

# Contrasts by castration group
  contrast_Genus_weaning= makeContrasts(weaningTruck-weaningFence_line,
                                       weaningFence_line-weaningNot_weaned,
                                       weaningTruck-weaningNot_weaned,
                                        levels=Mod_Genus_weaning)
  res_Genus_weaning = contrasts.fit(ZIGFit_Genus_weaning, contrast_Genus_weaning)
  resEB_Genus_weaning = eBayes(res_Genus_weaning)

# Extract results
  fz_Genus_weaning <- topTable(resEB_Genus_weaning, adjust.method="BH",number = 1200)
  
  fz_Genus_weaning_T_F <- topTable(resEB_Genus_weaning, adjust.method="BH",coef = 1,number = 1200)
  fz_Genus_weaning_F_Not <- topTable(resEB_Genus_weaning, adjust.method="BH",coef = 2,number = 1200)
  fz_Genus_weaning_T_Not <- topTable(resEB_Genus_weaning, adjust.method="BH",coef = 3,number = 1200)

# Write CSV files
  write.csv(fz_Genus_weaning_T_F, "./fz_Genus_weaning_T_F.csv")
  write.csv(fz_Genus_weaning_F_Not, "./fz_Genus_weaning_F_Not.csv")
  write.csv(fz_Genus_weaning_T_Not, "./fz_Genus_weaning_T_Not.csv")




```

## 7.2. PLOTTING Vulcano plots

```{r}

#1. WEANING PHYLUM

  log2.ALLDays<-read.csv(file = "fz_FINAL_weaning_FL-T.csv")  ##filename--contains model output
  # Subsetting phylum
  log2.phy.ALLDays <- subset(log2.ALLDays, log2.ALLDays$taxa_level == "gene_group")
  View(log2.ALLDays)
  
    #Subsetting weaning
  Weaning_log2.phy.ALLDays <- subset(log2.phy.ALLDays, log2.phy.ALLDays$variable == "weaning")
  Weaning_log2.phy.ALLDays <- filter(Weaning_log2.phy.ALLDays, taxa != "NA", taxa !="no_match")
  #median(log2.phy$AveExpr)
  Weaning_log2.phy.ALLDays$comparison_type=factor(Weaning_log2.phy.ALLDays$comparison_type, levels=c("Fence-line vs Not weaned", "Truck vs Not weaned", "Fence-line vs Truck"))

  log2.ALLDays$Significant <- ifelse(log2.ALLDays$adj.P.Val < 0.05, "AdjP < 0.05", "Not Sig")

# PLOT!
  vol_phy.ALLDays_weaning <- ggplot(log2.ALLDays, aes(x = logFC, y = taxa, fill = Significant))+
  #geom_point(aes(size=size_AveEx),width = 0.30, height = 0.35, alpha = 1, na.rm = T, shape = 21, colour = "black") +
  geom_point(aes(size=AveExpr), alpha = 1, na.rm = T, shape = 21, colour = "black") +
  scale_fill_manual(values = c("red", "grey70"))+                    #red significant 
  #theme(axis.ticks = element_line(color = "grey", size=(1)))
  #theme_hc()+
  theme(
    panel.grid.major = element_line(colour="#f0f0f0"),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=10, margin = margin(.5,0,.5,0)),
    strip.text.y=element_text(size=10, margin = margin(.5,0,.5,0)),
    axis.text.y=element_text(size=10),
    axis.title=element_text(size=10),
    legend.position="right",
    panel.spacing=unit(0.3, "lines"),
    plot.title=element_text(size=12, hjust=0.5),
    legend.text=element_text(size=8),
    legend.title=element_text(size=10),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA,color = "grey")) +
  xlim(-1.5,1.5)+
  geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(~collection_day, nrow=1)+
  ylab("")+
  ggtitle("Fence-line vs Truck weaning each collection day: Gene-level pairwise comparison")+         # change the title based on contrast 1,2..
  theme(legend.position = 'none')

vol_phy.ALLDays_weaning


    #Subsetting collection_day
  Day_log2.phy.ALLDays <- subset(log2.phy.ALLDays, log2.phy.ALLDays$variable == "collection_day")
  Day_log2.phy.ALLDays <- filter(Day_log2.phy.ALLDays, taxa != "NA", taxa !="no_match")
  #median(log2.phy$AveExpr)
  Day_log2.phy.ALLDays$comparison_type=factor(Day_log2.phy.ALLDays$comparison_type, levels=c("At weaning vs Pre-weaning", "Post-weaning vs Pre-weaning", "Post-weaning vs At weaning"))

  Day_log2.phy.ALLDays$Significant <- ifelse(Day_log2.phy.ALLDays$adj.P.Val < 0.05, "AdjP < 0.05", "Not Sig")

# PLOT!
  vol_phy.ALLDays_Day <- ggplot(Day_log2.phy.ALLDays, aes(x = logFC, y = taxa, fill = Significant))+
  #geom_point(aes(size=size_AveEx),width = 0.30, height = 0.35, alpha = 1, na.rm = T, shape = 21, colour = "black") +
  geom_point(aes(size=AveExpr), alpha = 1, na.rm = T, shape = 21, colour = "black") +
  scale_fill_manual(values = c("red", "grey70"))+                    #red significant 
  #theme(axis.ticks = element_line(color = "grey", size=(1)))
  #theme_hc()+
  theme(
    panel.grid.major = element_line(colour="#f0f0f0"),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=10, margin = margin(.5,0,.5,0)),
    strip.text.y=element_text(size=10, margin = margin(.5,0,.5,0)),
    axis.text.y=element_text(size=10),
    axis.title=element_text(size=10),
    legend.position="right",
    panel.spacing=unit(0.3, "lines"),
    plot.title=element_text(size=12, hjust=0.5),
    legend.text=element_text(size=8),
    legend.title=element_text(size=10),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA,color = "grey")) +
  xlim(-1.5,1.5)+
  geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(~comparison_type, nrow=1)+
  ylab("")+
  ggtitle("Weaning strategies: Phylum-level pairwise comparison")+         # change the title based on contrast 1,2..
  theme(legend.position = 'none')

vol_phy.ALLDays_Day





```




