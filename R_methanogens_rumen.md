---
title: "Methane_calves"
author: "G Diaz"
date: "2022-10-24"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 0. Loading packages & stablishing working directory

```{r libraries}

  setwd("~/Desktop/Gerardo/UMN/PhD/Methane/Prelim Data")
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

# Set sed number for random-based functions
  set.seed(143)

```

# 1. Loading raw data
Take into account that results come from a rarefaction process at minimun number of reads of the sample set: *6188129 reads* (check sampleinfo.txt)
Tax and otu tables have been manually built using the raw count table obtained from bioinformatic analysis ("function_profile" & "taxonomy_profile")

## 1.1. FUNCTION

```{r data}

  otu_mat<- read_excel("MCycDB_function_profile.xlsx", sheet = "otu_table")
  tax_mat<- read_excel("MCycDB_function_profile.xlsx", sheet = "taxa_table")
  samples_df<- read_excel("RumenDATA_phyloseq_Ver04.xlsx", sheet = "Samples")
  
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

# FORMATTING DATA
    otu_mat <- otu_mat %>%
    tibble::column_to_rownames("OTU") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("OTU")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("samples")

  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  
# TRANSFORMING TO PHYLOSEQ OBJECT
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)
  
  Meth1 <- phyloseq(OTU, TAX, samples)
  Meth1

# SUBSETTING DATA
  Meth1_Mock <- subset_samples(Meth1, sample_type =="Mock")
  Meth1_Sample <- subset_samples(Meth1, sample_type =="Sample")
  
  Meth1_Sample1 <- subset_samples(Meth1_Sample, collection_day =="Pre_weaning")
  Meth1_Sample2 <- subset_samples(Meth1_Sample, collection_day =="At_weaning")
  Meth1_Sample3 <- subset_samples(Meth1_Sample, collection_day =="Post_weaning")

# SUBSETTING FOR ALPHA DIVERSITY
  samples_2<-subset(samples_df, samples_df$sample_type=="Sample")
  samples_2=data.frame(samples_2)
  samples_D1<-subset(samples_2, samples_2$collection_day=="Pre_weaning")
  samples_D1=data.frame(samples_D1)
  samples_D2<-subset(samples_2, samples_2$collection_day=="At_weaning")
  samples_D2=data.frame(samples_D2)
  samples_D3<-subset(samples_2, samples_2$collection_day=="Post_weaning")
  samples_D3=data.frame(samples_D3)

# FORMATTING DATA FOR ALPHA DIVERSITY PLOTS
  samples_2$castration_group=factor(samples_2$castration_group, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning"))
  samples_2$castration=factor(samples_2$castration, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning", "Not_castrated"))
  samples_2$weaning_group=factor(samples_2$weaning_group, levels=c("Fence_line", "Truck"))
  samples_2$weaning=factor(samples_2$weaning, levels=c("Fence_line", "Truck", "Not_weaned"))
  samples_2$collection_day=factor(samples_2$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))

# FORMATTING DATA FOR BETA DIVERSITY PLOTS
# Changing from phyloseq object to MRexperiment (metagenomeSeq format) for CSS normalization
rumen.metaseq <- phyloseq_to_metagenomeSeq(Meth1_Sample) # change here
rumen.metaseq.norm<- cumNorm(rumen.metaseq, p=cumNormStat(rumen.metaseq))
CSS_rumen.metaseq <- MRcounts(rumen.metaseq.norm, norm = TRUE)
CSS_Meth1_Sample <- merge_phyloseq(otu_table(CSS_rumen.metaseq,
taxa_are_rows=T),sample_data(Meth1_Sample),tax_table(Meth1_Sample))
  
```

## 1.2. TAXONOMY_genus


```{r data}

  otu_mat2<- read_excel("MCycDB_taxonomy_profile.xlsx", sheet = "otu_table_genus")
  tax_mat2<- read_excel("MCycDB_taxonomy_profile.xlsx", sheet = "tax_table_genus")

  # FORMATTING DATA
    otu_mat2 <- otu_mat2 %>%
    tibble::column_to_rownames("OTU") 
  tax_mat2 <- tax_mat2 %>% 
    tibble::column_to_rownames("OTU")

  otu_mat2 <- as.matrix(otu_mat2)
  tax_mat2 <- as.matrix(tax_mat2)
  
# TRANSFORMING TO PHYLOSEQ OBJECT
  OTU2 = otu_table(otu_mat2, taxa_are_rows = TRUE)
  TAX2 = tax_table(tax_mat2)
  
  Meth2 <- phyloseq(OTU2, TAX2, samples)
  Meth2

# SUBSETTING DATA
  Meth2_Mock <- subset_samples(Meth2, sample_type =="Mock")
  Meth2_Sample <- subset_samples(Meth2, sample_type =="Sample")
  
  Meth2_Sample1 <- subset_samples(Meth2_Sample, collection_day =="Pre_weaning")
  Meth2_Sample2 <- subset_samples(Meth2_Sample, collection_day =="At_weaning")
  Meth2_Sample3 <- subset_samples(Meth2_Sample, collection_day =="Post_weaning")

# FORMATTING DATA FOR BETA DIVERSITY PLOTS
# Changing from phyloseq object to MRexperiment (metagenomeSeq format) for CSS normalization
rumen.metaseq <- phyloseq_to_metagenomeSeq(Meth2_Sample) # change here
rumen.metaseq.norm<- cumNorm(rumen.metaseq, p=cumNormStat(rumen.metaseq))
CSS_rumen.metaseq <- MRcounts(rumen.metaseq.norm, norm = TRUE)
CSS_Meth2_Sample <- merge_phyloseq(otu_table(CSS_rumen.metaseq,
taxa_are_rows=T),sample_data(Meth2_Sample),tax_table(Meth2_Sample))
  
```

### 1.2.1. Relative ab_genus

```{r}

CSS_Meth2_Sample

# 1. Converting to Percentage and transforming Genus tax_glom object to long-format table
phy_relative3 <- transform_sample_counts(CSS_Meth2_Sample, function(x) (x / sum(x))*100 )
phy_relative_long3 <- psmelt(phy_relative3)

# 2. Creating a new column (variable) for the  long-format table, called "mean_relative_abund" 
phy_relative_long3 <- phy_relative_long3 %>%
  group_by(Genus) %>%
  mutate(mean_relative_abund = mean(Abundance))

# 3. Formatting the "Genus" column variable as names and the "mean_relative_abund" column variable as continuous numbers
phy_relative_long3$Genus <- as.character(phy_relative_long3$Genus)
phy_relative_long3$mean_relative_abund <- as.numeric(phy_relative_long3$mean_relative_abund)

# 4. Giving the name "Others (< 1%)" to all the OTUs in the "Genus" column that have less than 1 in the "Abundace" column in the long-format table 
phy_relative_long3$Genus[phy_relative_long3$Abundance < 1] <- "Others (< 1%)" 

# 5. PLOTTING by collection day, castration and weaning group

phy_relative_long3 %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")



```

## 1.3. TAXONOMY_specie


```{r data}

  otu_mat3<- read_excel("MCycDB_taxonomy_profile.xlsx", sheet = "otu_table_species")
  tax_mat3<- read_excel("MCycDB_taxonomy_profile.xlsx", sheet = "tax_table_species")

  # FORMATTING DATA
    otu_mat3 <- otu_mat3 %>%
    tibble::column_to_rownames("OTU") 
  tax_mat3 <- tax_mat3 %>% 
    tibble::column_to_rownames("OTU")

  otu_mat3 <- as.matrix(otu_mat3)
  tax_mat3 <- as.matrix(tax_mat3)
  
# TRANSFORMING TO PHYLOSEQ OBJECT
  OTU3 = otu_table(otu_mat3, taxa_are_rows = TRUE)
  TAX3 = tax_table(tax_mat3)
  
  Meth3 <- phyloseq(OTU3, TAX3, samples)
  Meth3

# SUBSETTING DATA
  Meth3_Mock <- subset_samples(Meth3, sample_type =="Mock")
  Meth3_Sample <- subset_samples(Meth3, sample_type =="Sample")
  
  Meth3_Sample1 <- subset_samples(Meth3_Sample, collection_day =="Pre_weaning")
  Meth3_Sample2 <- subset_samples(Meth3_Sample, collection_day =="At_weaning")
  Meth3_Sample3 <- subset_samples(Meth3_Sample, collection_day =="Post_weaning")

# FORMATTING DATA FOR BETA DIVERSITY PLOTS
# Changing from phyloseq object to MRexperiment (metagenomeSeq format) for CSS normalization
rumen.metaseq <- phyloseq_to_metagenomeSeq(Meth3_Sample) # change here
rumen.metaseq.norm<- cumNorm(rumen.metaseq, p=cumNormStat(rumen.metaseq))
CSS_rumen.metaseq <- MRcounts(rumen.metaseq.norm, norm = TRUE)
CSS_Meth3_Sample <- merge_phyloseq(otu_table(CSS_rumen.metaseq,
taxa_are_rows=T),sample_data(Meth3_Sample),tax_table(Meth3_Sample))
```

# 2. ALPHA DIVERSITY

## 2.1. FUNCTION

```{r alpha-plot}

# Obtaining alpha diversity for every sample at OTU Level
  genus_alpha.div= estimate_richness(Meth1, measures=c("Observed", "InvSimpson", "Shannon"))
  write.csv(genus_alpha.div, "Alpha_Meth1_Gene.csv")

  genus_alpha.div2= evenness(Meth1, 'pielou')
  write.csv(genus_alpha.div2, "Alpha-Evenness_Meth1_Gene.csv")

# PLOTTING
# plotting richness, shannon and Pielou by WEANING and Collection day
   ggplot(samples_2, aes(collection_day,Kra0.1_Obs_Genus, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ geom_jitter(position=position_dodge(0.8),size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Richness")+xlab("Collection day")

  ggplot(samples_2, aes(collection_day,Kra0.1_Shan_Genus, color=weaning_group))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ geom_jitter(position=position_dodge(0.8), size=2)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Shannon Index")+xlab("Collection day")
   
  ggplot(samples_2, aes(collection_day,Kra0.1_Even_Genus, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ylim(0, 1)+ geom_jitter(position=position_dodge(0.8), size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Pielou's Evenness Index")+xlab("Collection day")
     
# STATISTICS
# Model and type-III ANOVA

  md= samples_2
  
  Model9_Shan = lmer(Kra0.1_Shan_Genus ~ weaning_group*collection_day + castration + age_days + (1|cow_ID), data=md, REML = T) 

  summary(Model9_Shan)
  confint(Model9_Shan)
  anova(Model9_Shan)
  emmeans(Model9_Shan,pairwise~weaning_group|collection_day)



```

## 2.1. TAXONOMY_genus

```{r}

# Obtaining alpha diversity for every sample at OTU Level
  genus_alpha.div= estimate_richness(Meth2, measures=c("Observed", "InvSimpson", "Shannon"))
  write.csv(genus_alpha.div, "Alpha_Meth2_Gene.csv")

  genus_alpha.div2= evenness(Meth2, 'pielou')
  write.csv(genus_alpha.div2, "Alpha-Evenness_Meth2_Gene.csv")

# PLOTTING
# plotting richness, shannon and Pielou by WEANING and Collection day
   ggplot(samples_2, aes(collection_day,Kra2_Obs_Genus, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ geom_jitter(position=position_dodge(0.8),size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Richness")+xlab("Collection day")

  ggplot(samples_2, aes(collection_day,Kra2_Shan_Genus, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ geom_jitter(position=position_dodge(0.8), size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Shannon Index")+xlab("Collection day")
   
  ggplot(samples_2, aes(collection_day,Kra2_Even_Genus, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ylim(0, 1)+ geom_jitter(position=position_dodge(0.8), size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Pielou's Evenness Index")+xlab("Collection day")
     
  
# STATISTICS
# Model and type-III ANOVA

  Model10_Shan = lmer(Kra2_Shan_Genus ~ weaning+collection_day + castration + age_days + seq_K0.1classified + (1|cow_ID), data=md, REML = T) 

  summary(Model10_Shan)
  confint(Model10_Shan)
  anova(Model10_Shan)
  emmeans(Model10_Shan,pairwise~weaning)

```

## 2.1. TAXONOMY_specie

```{r}

# Obtaining alpha diversity for every sample at OTU Level
  genus_alpha.div= estimate_richness(Meth3, measures=c("Observed", "InvSimpson", "Shannon"))
  write.csv(genus_alpha.div, "Alpha_Meth3_Gene.csv")

  genus_alpha.div2= evenness(Meth3, 'pielou')
  write.csv(genus_alpha.div2, "Alpha-Evenness_Meth3_Gene.csv")

  # PLOTTING
# plotting richness, shannon and Pielou by WEANING and Collection day
   ggplot(samples_2, aes(collection_day,Kra3_Obs_Genus, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ geom_jitter(position=position_dodge(0.8),size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Richness")+xlab("Collection day")

  ggplot(samples_2, aes(collection_day,Kra3_Shan_Genus, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ geom_jitter(position=position_dodge(0.8), size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Shannon Index")+xlab("Collection day")
   
  ggplot(samples_2, aes(collection_day,Kra3_Even_Genus, color=weaning))+
                 geom_boxplot(position=position_dodge(0.8),lwd=1) +
   scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+
   theme_classic()+ylim(0, 1)+ geom_jitter(position=position_dodge(0.8), size=3)+
   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
   ylab("Pielou's Evenness Index")+xlab("Collection day")
     
  
# STATISTICS
# Model and type-III ANOVA

  Model11_Shan = lmer(Kra3_Shan_Genus ~ weaning+collection_day + castration + age_days + seq_K0.1classified + (1|cow_ID), data=md, REML = T) 

  summary(Model11_Shan)
  confint(Model11_Shan)
  anova(Model11_Shan)
  emmeans(Model11_Shan,pairwise~weaning)


```

# 3. BETA DIVERSITY

## 3.1. FUNCTION

```{r beta-plots}

# 1. Get the ordination from Normalized and genus-level glom phyloseq object 
  genus.ord <- ordinate(CSS_Meth1_Sample, "NMDS", "bray")
  # stress: 0.1684933

# 2. Plotting the ordination per WEANING and collection day
# relabeling
levels(sample_data(CSS_Meth1_Sample)$weaning) <- c('Fence_line', 'Truck', 'Before_weaning')

p<-plot_ordination(CSS_Meth1_Sample, genus.ord, type="sample", color="weaning", shape="collection_day")+
  geom_point(size=4.5, alpha=0.8)+scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+scale_shape_manual(values = c(15,16,17))+  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+theme_classic()+stat_ellipse()


p$layers
p$layers <- p$layers[-1]
p

# STATISTICS 

# 1. Inputting data

metadata= data.frame(sample_data(CSS_Meth1_Sample))
dist = vegdist(t(otu_table(CSS_Meth1_Sample)), method ="bray" )

# 2. Models
set.seed(143)
  permanova5.1 = adonis2(dist ~ weaning + collection_day  ,strata=metadata$collection_day,  data=metadata)  
  permanova5.2 = adonis2(dist ~ castration + collection_day  ,strata=metadata$collection_day ,  data=metadata) 
  permanova5.3 = adonis2(dist ~ collection_day + weaning + castration,strata=metadata$weaning,  data=metadata) 

  #final model
permanova5.1.2 = adonis2(dist ~ weaning_group+collection_day+castration ,   data=metadata , by="margin")
# Rumen_dist_clean = removed R0048 because did not allow to get ordination for collection_day 1 and improved the figure but DID NOT CHANGE THE DIST MATRIX NOR STATISTICS 
  
```

## 3.2. TAXONOMY_genus

```{r beta-plots2}

# 1. Get the ordination from Normalized and genus-level glom phyloseq object 
  genus.ord <- ordinate(CSS_Meth2_Sample, "NMDS", "bray")
  # stress: 0.2099721

# 2. Plotting the ordination per WEANING and collection day

p2<-plot_ordination(CSS_Meth2_Sample, genus.ord, type="sample", color="weaning", shape="collection_day")+
  geom_point(size=4.5, alpha=0.8)+scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+scale_shape_manual(values = c(15,16,17))+  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+theme_classic()

p2$layers
p2$layers <- p2$layers[-1]
p2

# STATISTICS 

# 1. Inputting data

metadata= data.frame(sample_data(CSS_Meth2_Sample))
dist = vegdist(t(otu_table(CSS_Meth2_Sample)), method ="bray" )

# 2. Models

  permanova6.1 = adonis2(dist ~ weaning + collection_day  ,strata=metadata$collection_day,  data=metadata)  
  permanova6.2 = adonis2(dist ~ castration + collection_day  ,strata=metadata$collection_day ,  data=metadata) 
  permanova6.3 = adonis2(dist ~ collection_day + weaning + castration,strata=metadata$weaning,  data=metadata) 


```


## 3.3. TAXONOMY_species

```{r beta-plots3}

# 1. Get the ordination from Normalized and genus-level glom phyloseq object 
  genus.ord <- ordinate(CSS_Meth3_Sample, "PCoA", "bray")
  # stress: 8.284773e-05 

# 2. Plotting the ordination per WEANING and collection day

p3<-plot_ordination(CSS_Meth3_Sample, genus.ord, type="sample", color="weaning", shape="collection_day")+
  geom_point(size=4.5, alpha=0.8)+scale_fill_manual(values = c("#967bba", "#505359", "#c9c8c8"))+scale_shape_manual(values = c(15,16,17))+  scale_color_manual(values = c("#967bba", "#505359", "#c9c8c8"))+theme_classic()

p3$layers
p3$layers <- p3$layers[-1]
p3

# STATISTICS 

# 1. Inputting data

metadata= data.frame(sample_data(CSS_Meth3_Sample))
dist = vegdist(t(otu_table(CSS_Meth3_Sample)), method ="bray" )

# 2. Models

  permanova7.1 = adonis2(dist ~ weaning + collection_day  ,strata=metadata$collection_day,  data=metadata)  
  permanova7.2 = adonis2(dist ~ castration + collection_day  ,strata=metadata$collection_day ,  data=metadata) 
  permanova7.3 = adonis2(dist ~ collection_day + weaning + castration,strata=metadata$weaning,  data=metadata) 


```


#4. DIFFERENTIAL ABUNDANCE
## 4.1. FUNCTION

```{r diff}

# 1. DEFINE YOUR METAGENOMESEQ EXPERIMENT (cumNorm-ed)
# we need to filter the OTUs at Phylum level to be present in at least in 2 samples (see ?filterData for help)
  Rumen_metaseq <- phyloseq_to_metagenomeSeq(CSS_Meth1_Sample) # it will change to MRexperiment as required by the "filterData" function
  Rumen_Gene_filtered <- filterData(aggTax(Rumen_metaseq, lvl="Gene"), present = 2) # it will aggregate taxa at phylum (and give phylum names instead of only "OTU XXXX"), and filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  Rumen_Gene_filtered.metaseq <- cumNorm(Rumen_Gene_filtered) # it will normalize again the new filtered MRexperiment 

# 2. DEFINE YOUR METADATA VARIABLES TO INCLUDE IN THE MODEL (pData)
  
  collection_day <- pData(Rumen_Gene_filtered.metaseq)$collection_day
  castration <-  pData(Rumen_Gene_filtered.metaseq)$castration # will include 4 castration strategies + Not_castrated
  weaning <-  pData(Rumen_Gene_filtered.metaseq)$weaning # will include 2 weaning strategies + Not_weaned
  cow_ID <- pData(Rumen_Gene_filtered.metaseq)$cow_ID

# 3. Define your settings (zigControl), we will use the ones by default: Same as before
  settings = zigControl(maxit = 10, verbose = TRUE)

# 4. DEFINE YOUR MODEL ZERO TO ACCOUNT FOR DEPTH OF COVERAGE (zero_mod)
## MODEL 0: ZERO MODEL FOR ALL DAYS (depth of coverage of the NOT-normalized (RAW) MRexperiment)
  
  rumen.metaseq <- phyloseq_to_metagenomeSeq(Meth1_Sample)
  Rumen_RAW_Gene_filtered <- filterData(aggTax(rumen.metaseq,lvl="Gene"), present = 2) # it will filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  zero_mod_Gene <- model.matrix(~0+log(libSize(Rumen_RAW_Gene_filtered))) # model it!
  
# 5. DEFINE YOUR MODEL TO TEST (model.matrix), FIT THE MODEL (fitZig), MAKE PAIRWISE CONTRASTS (makeContrasts, contrasts.fit, eBayes) & EXTRACT RESULTS (topTable) 

# MODEL 2: Weaning overall

# Assign model (model.matrix)
  design_weaning = model.matrix(~0 + weaning)
  #dup <- duplicateCorrelation(MRcounts(Rumen_Gene_filtered.metaseq),block=cow_ID, design= design_weaning)
  
# Assign correlated data (random effect)
  v_weaning <- voom(MRcounts(Rumen_Gene_filtered.metaseq), design_weaning)
  dup_weaning <- duplicateCorrelation(v_weaning, block=cow_ID, design= design_weaning)
  dup_weaning$consensus.correlation # to check your correlation coeff?

# Make ZIG model 
################## START PIECE OF CODE DIDN'T WORK #########################
# DOES NOT WORK BECAUSE IT IS LIMMA VERSION 3.52.4 FOR R 4.2.1
# IT WORKED IN LIMMA VERSION 3.50.3 FOR R 4.1.0 (Windows)
  ZIG_Gene_weaning = fitZig(obj= Rumen_Gene_filtered.metaseq, 
                              mod = design_weaning, 
                              control = settings, 
                              zeroMod=zero_mod_Gene, 
                              useCSSoffset=FALSE,
                              useMixedModel=dup_weaning$consensus.correlation, 
                              block=cow_ID)
  
  ZIGFit_Gene_weaning = ZIG_Gene_weaning@fit
  Mod_Gene_weaning = ZIG_Gene_weaning@fit$design

# Contrasts by weaning group
  contrast_Gene_weaning= makeContrasts(weaningFence_line-weaningTruck,
                                       weaningNot_weaned-weaningFence_line,
                                       weaningNot_weaned-weaningTruck,
                                        levels=Mod_Gene_weaning)
  res_Gene_weaning = contrasts.fit(ZIGFit_Gene_weaning, contrast_Gene_weaning)
  resEB_Gene_weaning = eBayes(res_Gene_weaning)

# Extract results
  fz_Gene_weaning <- topTable(resEB_Gene_weaning, adjust.method="BH",number = 1200)
  
  fz_Gene_weaning_F_T <- topTable(resEB_Gene_weaning, adjust.method="BH",coef = 1,number = 12000)
  fz_Gene_weaning_Not_F <- topTable(resEB_Gene_weaning, adjust.method="BH",coef = 2,number = 12000)
  fz_Gene_weaning_Not_T <- topTable(resEB_Gene_weaning, adjust.method="BH",coef = 3,number = 12000)

# Write CSV files
  write.csv(fz_Gene_weaning_F_T, "./fz_Gene_weaning_F_T.csv")
  write.csv(fz_Gene_weaning_Not_F, "./fz_Gene_weaning_Not_F.csv")
  write.csv(fz_Gene_weaning_Not_T, "./fz_Gene_weaning_Not_T.csv")

####### FINISHED PIECE OF CODE DIDN'T WORK ################# 
  

# VULCANO PLOT
log2.ALLDays<-read.csv(file = "LogFold_Rumen_ALLDays.csv")  ##filename--contains model output
# Subsetting gene
  log2.gene.ALLDays <- subset(log2.ALLDays, log2.ALLDays$taxa_level == "gene")
  log2.gene.ALLDays$Significant <- ifelse(log2.gene.ALLDays$adj.P.Val < 0.05, "AdjP < 0.05", "Not Sig")
#log2.phy.Day3$EffectSize <- ifelse(log2.phy.Day3$AveExpr < 4.2057385, "ES<1", "EF>1")
#topTable_contrasts_df['taxa'] <- row.names(topTable_contrasts_df)
#log2.phy.Day3_fil <- subset(log2.phy.Day3,  percentile_fifty_cut%in% c("1"))       ## taxa av.exp >50th percentile
#log2.phy.Day3_fil

#size_AveEx <- sqrt(log2.phy.Day3_fil$AveExpr/pi)    ## size of average expression/effectsize

  View(log2.gene.ALLDays)
  
  #Subsetting weaning comparisons
  Weaning_log2.gene.ALLDays_Not_F <- subset(log2.gene.ALLDays, log2.gene.ALLDays$comparison_type == "Not weaned vs Fence-line")
  Weaning_log2.gene.ALLDays_Not_T <- subset(log2.gene.ALLDays, log2.gene.ALLDays$comparison_type == "Not weaned vs Truck")
  Weaning_log2.gene.ALLDays_F_T <- subset(log2.gene.ALLDays, log2.gene.ALLDays$comparison_type == "Fence-line vs Truck")
  #median(log2.phy$AveExpr)

# PLOT: Not weaned vs Fence-line
  
  Not_F3 <- ggplot(Weaning_log2.gene.ALLDays_Not_F, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.gene.ALLDays_Not_F,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))
  
# PLOT: Not weaned vs Truck
  
  Not_T3 <- ggplot(Weaning_log2.gene.ALLDays_Not_T, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.gene.ALLDays_Not_T,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))

# PLOT: Fence-line vs Truck
# 54 genes increased in Fence-line (compared to truck)
# 108 gene decreased in Fence-line (compared to truck)
  
  F_T3 <- ggplot(Weaning_log2.gene.ALLDays_F_T, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.gene.ALLDays_F_T,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))



############# DELETE #######

#Subsetting weaning
Weaning_log2.gene.ALLDays <- subset(log2.gene.ALLDays, log2.gene.ALLDays$variable == "weaning")
Weaning_log2.gene.ALLDays <- filter(Weaning_log2.gene.ALLDays, taxa != "NA", taxa !="no_match")
#median(log2.phy$AveExpr)
Weaning_log2.gene.ALLDays$comparison_type=factor(Weaning_log2.gene.ALLDays$comparison_type, levels=c("Not weaned vs Fence-line", "Not weaned vs Truck", "Fence-line vs Truck"))

Weaning_log2.gene.ALLDays$Significant <- ifelse(Weaning_log2.gene.ALLDays$adj.P.Val < 0.05, "AdjP < 0.05", "Not Sig")

View(log2.gen.ALLDays)




```


## 4.2. TAXONOMY

```{r diff2}

# Phylum level 
# 0. Call and normalize phyloseq object

  otu_mat4<- read_excel("MCycDB_taxonomy_profile.xlsx", sheet = "otu_table_phylum")
  tax_mat4<- read_excel("MCycDB_taxonomy_profile.xlsx", sheet = "tax_table_phylum")

  # FORMATTING DATA
    otu_mat4 <- otu_mat4 %>%
    tibble::column_to_rownames("OTU") 
  tax_mat4 <- tax_mat4 %>% 
    tibble::column_to_rownames("OTU")

  otu_mat4 <- as.matrix(otu_mat4)
  tax_mat4 <- as.matrix(tax_mat4)
  
# TRANSFORMING TO PHYLOSEQ OBJECT
  OTU4 = otu_table(otu_mat4, taxa_are_rows = TRUE)
  TAX4 = tax_table(tax_mat4)
  
  Meth4 <- phyloseq(OTU4, TAX4, samples)
  Meth4

# SUBSETTING DATA
  Meth4_Sample <- subset_samples(Meth4, sample_type =="Sample")

# FORMATTING DATA FOR BETA DIVERSITY PLOTS
# Changing from phyloseq object to MRexperiment (metagenomeSeq format) for CSS normalization
rumen.metaseq <- phyloseq_to_metagenomeSeq(Meth4_Sample) # change here
rumen.metaseq.norm<- cumNorm(rumen.metaseq, p=cumNormStat(rumen.metaseq))
CSS_rumen.metaseq <- MRcounts(rumen.metaseq.norm, norm = TRUE)
CSS_Meth4_Sample <- merge_phyloseq(otu_table(CSS_rumen.metaseq,
taxa_are_rows=T),sample_data(Meth4_Sample),tax_table(Meth4_Sample))
  

# 1. DEFINE YOUR METAGENOMESEQ EXPERIMENT (cumNorm-ed)
# we need to filter the OTUs at Phylum level to be present in at least in 2 samples (see ?filterData for help)
  Rumen_metaseq <- phyloseq_to_metagenomeSeq(CSS_Meth4_Sample) # it will change to MRexperiment as required by the "filterData" function
  Rumen_Phy_filtered <- filterData(aggTax(Rumen_metaseq, lvl="Phylum"), present = 2) # it will aggregate taxa at phylum (and give phylum names instead of only "OTU XXXX"), and filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  Rumen_Phy_filtered.metaseq <- cumNorm(Rumen_Phy_filtered) # it will normalize again the new filtered MRexperiment 

# 2. DEFINE YOUR METADATA VARIABLES TO INCLUDE IN THE MODEL (pData)
  
  collection_day <- pData(Rumen_Phy_filtered.metaseq)$collection_day
  castration <-  pData(Rumen_Phy_filtered.metaseq)$castration # will include 4 castration strategies + Not_castrated
  weaning <-  pData(Rumen_Phy_filtered.metaseq)$weaning # will include 2 weaning strategies + Not_weaned
  cow_ID <- pData(Rumen_Phy_filtered.metaseq)$cow_ID

# 3. Define your settings (zigControl), we will use the ones by default: Same as before
  settings = zigControl(maxit = 10, verbose = TRUE)

# 4. DEFINE YOUR MODEL ZERO TO ACCOUNT FOR DEPTH OF COVERAGE (zero_mod)
## MODEL 0: ZERO MODEL FOR ALL DAYS (depth of coverage of the NOT-normalized (RAW) MRexperiment)
  
  rumen.metaseq # We will use this MRexperiment (metagenomeseq object) from previous normalization process. This file contains 95 samples for 3 days and all rank levels
  Rumen_RAW_Phylum_filtered <- filterData(aggTax(rumen.metaseq,lvl="Phylum"), present = 2) # it will filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  zero_mod_Phy <- model.matrix(~0+log(libSize(Rumen_RAW_Phylum_filtered))) # model it!
  
# 5. DEFINE YOUR MODEL TO TEST (model.matrix), FIT THE MODEL (fitZig), MAKE PAIRWISE CONTRASTS (makeContrasts, contrasts.fit, eBayes) & EXTRACT RESULTS (topTable) 

# MODEL 2: Weaning overall

# Assign model (model.matrix)
  design_weaning = model.matrix(~0 + weaning)
  #dup <- duplicateCorrelation(MRcounts(Rumen_Phy_filtered.metaseq),block=cow_ID, design= design_weaning)
  
# Assign correlated data (random effect)
  v_weaning <- voom(MRcounts(Rumen_Phy_filtered.metaseq), design_weaning)
  dup_weaning <- duplicateCorrelation(v_weaning, block=cow_ID, design= design_weaning)
  dup_weaning$consensus.correlation # to check your correlation coeff?

# Make ZIG model 
################## START PIECE OF CODE DIDN'T WORK #########################
# DOES NOT WORK BECAUSE IT IS LIMMA VERSION 3.52.4 FOR R 4.2.1
# IT WORKED IN LIMMA VERSION 3.50.3 FOR R 4.1.0 (Windows)
  ZIG_Phylum_weaning = fitZig(obj= Rumen_Phy_filtered.metaseq, 
                              mod = design_weaning, 
                              control = settings, 
                              zeroMod=zero_mod_Phy, 
                              useCSSoffset=FALSE,
                              useMixedModel=dup_weaning$consensus.correlation, 
                              block=cow_ID)
  
  ZIGFit_Phylum_weaning = ZIG_Phylum_weaning@fit
  Mod_Phylum_weaning = ZIG_Phylum_weaning@fit$design

# Contrasts by weaning group
  contrast_Phylum_weaning= makeContrasts(weaningFence_line-weaningTruck,
                                       weaningNot_weaned-weaningFence_line,
                                       weaningNot_weaned-weaningTruck,
                                        levels=Mod_Phylum_weaning)
  res_Phylum_weaning = contrasts.fit(ZIGFit_Phylum_weaning, contrast_Phylum_weaning)
  resEB_Phylum_weaning = eBayes(res_Phylum_weaning)

# Extract results
  fz_Phylum_weaning <- topTable(resEB_Phylum_weaning, adjust.method="BH",number = 1200)
  
  fz_Phylum_weaning_F_T <- topTable(resEB_Phylum_weaning, adjust.method="BH",coef = 1,number = 12000)
  fz_Phylum_weaning_Not_F <- topTable(resEB_Phylum_weaning, adjust.method="BH",coef = 2,number = 12000)
  fz_Phylum_weaning_Not_T <- topTable(resEB_Phylum_weaning, adjust.method="BH",coef = 3,number = 12000)

# Write CSV files
  write.csv(fz_Phylum_weaning_F_T, "./fz_Phylum_weaning_F_T.csv")
  write.csv(fz_Phylum_weaning_Not_F, "./fz_Phylum_weaning_Not_F.csv")
  write.csv(fz_Phylum_weaning_Not_T, "./fz_Phylum_weaning_Not_T.csv")

####### FINISHED PIECE OF CODE DIDN'T WORK ################# 
  

# Genus level

# 1. DEFINE YOUR METAGENOMESEQ EXPERIMENT (cumNorm-ed)
# we need to filter the OTUs at Genus level to be present in at least in 2 samples (see ?filterData for help). We will use previous MRexperiment used to filter at the Genus level:
  Rumen_metaseq <- phyloseq_to_metagenomeSeq(CSS_Meth2_Sample) # it will change to MRexperiment as required by the "filterData" function
  Rumen_Gen_filtered <- filterData(aggTax(Rumen_metaseq, lvl="Genus"), present = 2) # it will aggregate taxa at Genus (and give Genus names instead of only "OTU XXXX"), and filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  Rumen_Gen_filtered.metaseq <- cumNorm(Rumen_Gen_filtered) # it will normalize again the new filtered MRexperiment 

# 2. DEFINE YOUR MODEL ZERO TO ACCOUNT FOR DEPTH OF COVERAGE (zero_mod)
## MODEL 0: ZERO MODEL FOR ALL DAYS (depth of coverage of the NOT-normalized OTU table from glom phyloseq obj)
rumen.metaseq <- phyloseq_to_metagenomeSeq(Meth2_Sample) # change here
  Rumen_RAW_Genus_filtered <- filterData(aggTax(rumen.metaseq,lvl="Genus"), present = 2) # it will filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  zero_mod_Gen <- model.matrix(~0+log(libSize(Rumen_RAW_Genus_filtered))) # model it!

# MODEL 3: Weaning overall

# Assign model (model.matrix): same as Phylum
  design_weaning
# Assign correlated data (random effect)
  v_weaning_gen <- voom(MRcounts(Rumen_Gen_filtered.metaseq), design_weaning)

  dup_weaning_gen <- duplicateCorrelation(v_weaning_gen, block=cow_ID, design= design_weaning)

  dup_weaning_gen$consensus.correlation # to check your correlation coeff?

# Make ZIG model 
################## START PIECE OF CODE DIDN'T WORK #########################
# DOES NOT WORK BECAUSE IT IS LIMMA VERSION 3.52.4 FOR R 4.2.1
# IT WORKED IN LIMMA VERSION 3.50.3 FOR R 4.1.0 (Windows)
  ZIG_Genus_weaning = fitZig(obj= Rumen_Gen_filtered.metaseq, 
                             mod = design_weaning, 
                             control = settings, 
                             zeroMod=zero_mod_Gen, useCSSoffset=FALSE, 
                             useMixedModel=dup_weaning_gen$consensus.correlation, 
                             block=cow_ID)
  ZIGFit_Genus_weaning = ZIG_Genus_weaning@fit
  Mod_Genus_weaning = ZIG_Genus_weaning@fit$design

# Contrasts by weaning group
  contrast_Genus_weaning= makeContrasts(weaningFence_line-weaningTruck,
                                       weaningNot_weaned-weaningFence_line,
                                       weaningNot_weaned-weaningTruck,
                                        levels=Mod_Genus_weaning)
  res_Genus_weaning = contrasts.fit(ZIGFit_Genus_weaning, contrast_Genus_weaning)
  resEB_Genus_weaning = eBayes(res_Genus_weaning)

# Extract results
  fz_Genus_weaning <- topTable(resEB_Genus_weaning, adjust.method="BH",number = 12000)
  
  fz_Genus_weaning_F_T <- topTable(resEB_Genus_weaning, adjust.method="BH",coef = 1,number = 12000)
  fz_Genus_weaning_Not_F <- topTable(resEB_Genus_weaning, adjust.method="BH",coef = 2,number = 12000)
  fz_Genus_weaning_Not_T <- topTable(resEB_Genus_weaning, adjust.method="BH",coef = 3,number = 12000)

# Write CSV files
  write.csv(fz_Genus_weaning_F_T, "./fz_Genus_weaning_F_T.csv")
  write.csv(fz_Genus_weaning_Not_F, "./fz_Genus_weaning_Not_F.csv")
  write.csv(fz_Genus_weaning_Not_T, "./fz_Genus_weaning_Not_T.csv")

####### FINISHED PIECE OF CODE DIDN'T WORK ################# 
  
  
  
# SPECIES LEVEL
# 1. DEFINE YOUR METAGENOMESEQ EXPERIMENT (cumNorm-ed)
# we need to filter the OTUs at Genus level to be present in at least in 2 samples (see ?filterData for help). We will use previous MRexperiment used to filter at the Genus level:
  Rumen_metaseq <- phyloseq_to_metagenomeSeq(CSS_Meth3_Sample) # it will change to MRexperiment as required by the "filterData" function
  Rumen_Spe_filtered <- filterData(aggTax(Rumen_metaseq, lvl="Species"), present = 2) # it will aggregate taxa at Genus (and give Genus names instead of only "OTU XXXX"), and filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  Rumen_Spe_filtered.metaseq <- cumNorm(Rumen_Spe_filtered) # it will normalize again the new filtered MRexperiment 

# 2. DEFINE YOUR MODEL ZERO TO ACCOUNT FOR DEPTH OF COVERAGE (zero_mod)
## MODEL 0: ZERO MODEL FOR ALL DAYS (depth of coverage of the NOT-normalized OTU table from glom phyloseq obj)
rumen.metaseq <- phyloseq_to_metagenomeSeq(Meth3_Sample) # change here
  Rumen_RAW_Species_filtered <- filterData(aggTax(rumen.metaseq,lvl="Species"), present = 2) # it will filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  zero_mod_Spe <- model.matrix(~0+log(libSize(Rumen_RAW_Species_filtered))) # model it!

# MODEL 3: Weaning overall

# Assign model (model.matrix): same as Phylum
  design_weaning
# Assign correlated data (random effect)
  v_weaning_spe <- voom(MRcounts(Rumen_Spe_filtered.metaseq), design_weaning)

  dup_weaning_spe <- duplicateCorrelation(v_weaning_spe, block=cow_ID, design= design_weaning)

  dup_weaning_spe$consensus.correlation # to check your correlation coeff?

# Make ZIG model 
################## START PIECE OF CODE DIDN'T WORK #########################
# DOES NOT WORK BECAUSE IT IS LIMMA VERSION 3.52.4 FOR R 4.2.1
# IT WORKED IN LIMMA VERSION 3.50.3 FOR R 4.1.0 (Windows)
  ZIG_Species_weaning = fitZig(obj= Rumen_Spe_filtered.metaseq, 
                             mod = design_weaning, 
                             control = settings, 
                             zeroMod=zero_mod_Spe, useCSSoffset=FALSE, 
                             useMixedModel=dup_weaning_spe$consensus.correlation, 
                             block=cow_ID)
  ZIGFit_Species_weaning = ZIG_Species_weaning@fit
  Mod_Species_weaning = ZIG_Species_weaning@fit$design

# Contrasts by weaning group
  contrast_Species_weaning= makeContrasts(weaningFence_line-weaningTruck,
                                       weaningNot_weaned-weaningFence_line,
                                       weaningNot_weaned-weaningTruck,
                                        levels=Mod_Species_weaning)
  res_Species_weaning = contrasts.fit(ZIGFit_Species_weaning, contrast_Species_weaning)
  resEB_Species_weaning = eBayes(res_Species_weaning)

# Extract results
  fz_Species_weaning <- topTable(resEB_Species_weaning, adjust.method="BH",number = 12000)
  
  fz_Species_weaning_F_T <- topTable(resEB_Species_weaning, adjust.method="BH",coef = 1,number = 12000)
  fz_Species_weaning_Not_F <- topTable(resEB_Species_weaning, adjust.method="BH",coef = 2,number = 12000)
  fz_Species_weaning_Not_T <- topTable(resEB_Species_weaning, adjust.method="BH",coef = 3,number = 12000)

# Write CSV files
  write.csv(fz_Species_weaning_F_T, "./fz_Species_weaning_F_T.csv")
  write.csv(fz_Species_weaning_Not_F, "./fz_Species_weaning_Not_F.csv")
  write.csv(fz_Species_weaning_Not_T, "./fz_Species_weaning_Not_T.csv")

  
################## VULCANO PLOTS #################
  

#1. WEANING PHYLUM

  log2.ALLDays<-read.csv(file = "LogFold_Rumen_ALLDays.csv")  ##filename--contains model output
  # Subsetting phylum
  log2.phy.ALLDays <- subset(log2.ALLDays, log2.ALLDays$taxa_level == "phylum")
  View(log2.phy.ALLDays)
  
    #Subsetting weaning
  Weaning_log2.phy.ALLDays <- subset(log2.phy.ALLDays, log2.phy.ALLDays$variable == "weaning")
  Weaning_log2.phy.ALLDays <- filter(Weaning_log2.phy.ALLDays, taxa != "NA", taxa !="no_match")
  #median(log2.phy$AveExpr)
  Weaning_log2.phy.ALLDays$comparison_type=factor(Weaning_log2.phy.ALLDays$comparison_type, levels=c("Not weaned vs Fence-line", "Not weaned vs Truck", "Fence-line vs Truck"))

  Weaning_log2.phy.ALLDays$Significant <- ifelse(Weaning_log2.phy.ALLDays$adj.P.Val < 0.05, "AdjP < 0.05", "Not Sig")

# PLOT!
  vol_phy.ALLDays_weaning <- ggplot(Weaning_log2.phy.ALLDays, aes(x = logFC, y = taxa, fill = Significant))+
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

vol_phy.ALLDays_weaning
  

#2. WEANING GENUS

# Subsetting genus
  log2.gen.ALLDays <- subset(log2.ALLDays, log2.ALLDays$taxa_level == "genus")
  log2.gen.ALLDays$Significant <- ifelse(log2.gen.ALLDays$adj.P.Val < 0.05, "AdjP < 0.05", "Not Sig")
#log2.phy.Day3$EffectSize <- ifelse(log2.phy.Day3$AveExpr < 4.2057385, "ES<1", "EF>1")
#topTable_contrasts_df['taxa'] <- row.names(topTable_contrasts_df)
#log2.phy.Day3_fil <- subset(log2.phy.Day3,  percentile_fifty_cut%in% c("1"))       ## taxa av.exp >50th percentile
#log2.phy.Day3_fil

#size_AveEx <- sqrt(log2.phy.Day3_fil$AveExpr/pi)    ## size of average expression/effectsize

  View(log2.gen.ALLDays)
  
  #Subsetting weaning comparisons
  Weaning_log2.gen.ALLDays_Not_F <- subset(log2.gen.ALLDays, log2.gen.ALLDays$comparison_type == "Not weaned vs Fence-line")
  Weaning_log2.gen.ALLDays_Not_T <- subset(log2.gen.ALLDays, log2.gen.ALLDays$comparison_type == "Not weaned vs Truck")
  Weaning_log2.gen.ALLDays_F_T <- subset(log2.gen.ALLDays, log2.gen.ALLDays$comparison_type == "Fence-line vs Truck")
  #median(log2.phy$AveExpr)

# PLOT: Not weaned vs Fence-line
  
  Not_F <- ggplot(Weaning_log2.gen.ALLDays_Not_F, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.gen.ALLDays_Not_F,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))
  
# PLOT: Not weaned vs Truck
  
  Not_T <- ggplot(Weaning_log2.gen.ALLDays_Not_T, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.gen.ALLDays_Not_T,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))

# PLOT: Fence-line vs Truck
# 45 genera increased in Fence-line (compared to truck)
# 64 genera decreased in Fence-line (compared to truck)
  
  F_T <- ggplot(Weaning_log2.gen.ALLDays_F_T, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.gen.ALLDays_F_T,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))+
  ggtitle("Post-weaning: Phylum-level pairwise comparison")        # change the title based on contrast 1,2..


#2. WEANING SPECIES

# Subsetting genus
  log2.spe.ALLDays <- subset(log2.ALLDays, log2.ALLDays$taxa_level == "species")
  log2.spe.ALLDays$Significant <- ifelse(log2.spe.ALLDays$adj.P.Val < 0.05, "AdjP < 0.05", "Not Sig")
#log2.phy.Day3$EffectSize <- ifelse(log2.phy.Day3$AveExpr < 4.2057385, "ES<1", "EF>1")
#topTable_contrasts_df['taxa'] <- row.names(topTable_contrasts_df)
#log2.phy.Day3_fil <- subset(log2.phy.Day3,  percentile_fifty_cut%in% c("1"))       ## taxa av.exp >50th percentile
#log2.phy.Day3_fil

#size_AveEx <- sqrt(log2.phy.Day3_fil$AveExpr/pi)    ## size of average expression/effectsize

  View(log2.spe.ALLDays)
  
  #Subsetting weaning comparisons
  Weaning_log2.spe.ALLDays_Not_F <- subset(log2.spe.ALLDays, log2.spe.ALLDays$comparison_type == "Not weaned vs Fence-line")
  Weaning_log2.spe.ALLDays_Not_T <- subset(log2.spe.ALLDays, log2.spe.ALLDays$comparison_type == "Not weaned vs Truck")
  Weaning_log2.spe.ALLDays_F_T <- subset(log2.spe.ALLDays, log2.spe.ALLDays$comparison_type == "Fence-line vs Truck")
  #median(log2.phy$AveExpr)

# PLOT: Not weaned vs Fence-line
  
  Not_F2 <- ggplot(Weaning_log2.spe.ALLDays_Not_F, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.spe.ALLDays_Not_F,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))
  
# PLOT: Not weaned vs Truck
  
  Not_T2 <- ggplot(Weaning_log2.spe.ALLDays_Not_T, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.spe.ALLDays_Not_T,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))

# PLOT: Truck vs Fence-line
  
  F_T2 <- ggplot(Weaning_log2.spe.ALLDays_F_T, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~ group_comp+treatment_group)+
  geom_text_repel(data=subset(Weaning_log2.spe.ALLDays_F_T,adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))


```





