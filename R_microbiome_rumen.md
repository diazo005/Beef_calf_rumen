---
title: "Rumen_kraken_ver06"
author: "Gerardo Diaz"
date: "`r Sys.Date()`"
output: html_document
editor_options: 
  chunk_output_type: console
---

# 0. Loading packages
Code to load all the needed packages:
```{r libraries, warning=FALSE, message=FALSE}
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
Define you 3 tables to build your phyloseq object:
* OTU
* Taxonomy
* Samples
```{r}
  otu_mat<- read_excel("RumenDATA_phyloseq_Ver07.xlsx", sheet = "OTUS_matrix0.1")
  tax_mat<- read_excel("RumenDATA_phyloseq_Ver07.xlsx", sheet = "Taxonomy_table0.1")
  samples_df<- read_excel("RumenDATA_phyloseq_Ver07.xlsx", sheet = "Samples")
```

Give some format to the sample table (re-leveling, etc)
```{r}
# formatting and re-leveling categorical variables 
  samples_df$collection_date= factor(samples_df$collection_date, levels=c("9.21.21", "10.18.21", "10.20.21"))
  samples_df$collection_day= factor(samples_df$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))
   samples_df$castration_group= factor(samples_df$castration_group, levels=c("Birth", "Turnout", "Pre_weaning","Weaning"))
   samples_df$castration= factor(samples_df$castration, levels=c("Birth", "Turnout", "Pre_weaning","Weaning", "Not_castrated"))
   samples_df$castration_status= factor(samples_df$castration_status, levels = c("Castrated", "Not_castrated"))
   samples_df$weaning= factor(samples_df$weaning, levels=c("Fence_line", "Truck","Before_weaning"))
   samples_df$weaning_group=factor(samples_df$weaning_group, levels=c("Fence_line", "Truck"))
   samples_df$weaning_status=factor(samples_df$weaning_status, levels=c("Weaned", "Not_weaned"))
  samples_df$cow_ID= as.factor(samples_df$cow_ID)
  samples_df$extraction_date= as.factor(samples_df$extraction_date)
  samples_df$extraction_run= as.factor(samples_df$extraction_run)
```

## 1.1. Phyloseq object
Create the RAW phyloseq object
```{r}
# Define the row names from the otu column in otu_mat and tax_mat
  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("OTU") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("OTU")
   samples_df <- samples_df %>% 
    tibble::column_to_rownames("samples")
# Transform into matrixes otu and tax tables (sample table can be left as data frame)
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
# Transform to phyloseq objects
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples_df)
  ps.object <- phyloseq(OTU, TAX, samples)
```

Give a name to the ps.object
```{r}
# Change the left part of the equivalence 
  Rumen = ps.object
```

* Rumen (kraken confidence 0.1): 10894 taxa

# 2. Summary statistics

Just to explore what is going on
```{r}
samples_2<-subset(samples_df, samples_df$sample_type=="Sample")
samples_2=data.frame(samples_2)
samples_D1<-subset(samples_2, samples_2$collection_day=="Pre_weaning")
samples_D1=data.frame(samples_D1)
samples_D2<-subset(samples_2, samples_2$collection_day=="At_weaning")
samples_D2=data.frame(samples_D2)
samples_D3<-subset(samples_2, samples_2$collection_day=="Post_weaning")
samples_D3=data.frame(samples_D3)

##################################################################
# kruskal.test(samples_D3$age_days~samples_D3$weaning)           #
# wilcox.test(samples_D3$age_days~samples_D3$weaning)            #
# t.test(samples_D3$age_days~samples_D3$weaning)                 #
# ANOVA<-aov(samples_D3$age_days~samples_D3$weaning)             #
# summary(ANOVA)                                                 #
##################################################################

# CASTRATION GROUP
describeBy(samples_D3$age_days,samples_D3$castration_group)
kruskal.test(samples_D3$age_days~samples_D3$castration_group)

describeBy(samples_D3$weight_birth,samples_D3$castration_group)
kruskal.test(samples_D3$weight_birth~samples_D3$castration_group)

describeBy(samples_D3$weight,samples_D3$castration_group)
kruskal.test(samples_D3$weight~samples_D3$castration_group)

describeBy(samples_D3$weight_ADGUpToDate,samples_D3$castration_group)
kruskal.test(samples_D3$weight_ADGUpToDate~samples_D3$castration_group)

# WEANING GROUP
describeBy(samples_D3$age_days,samples_D3$weaning_group)
kruskal.test(samples_D3$age_days~samples_D3$weaning_group)

describeBy(samples_D3$weight_birth,samples_D3$weaning_group)
kruskal.test(samples_D3$weight_birth~samples_D3$weaning_group)

describeBy(samples_D3$weight,samples_D3$weaning_group)
kruskal.test(samples_D3$weight~samples_D3$weaning_group)

describeBy(samples_D3$weight_ADGUpToDate,samples_D3$weaning_group)
kruskal.test(samples_D3$weight_ADGUpToDate~samples_D3$weaning_group)

# RAW SEQUENCING READS
describeBy(samples_2$seq_RAWreads, samples_2$collection_day)
kruskal.test(samples_2$seq_RAWreads, samples_2$collection_day)

# NON-HOST READS
describeBy(samples_2$seq_nonHOSTreads, samples_2$collection_day)
kruskal.test(samples_2$seq_nonHOSTreads, samples_2$collection_day)

# CLASSIFIED READS
describeBy(samples_df$seq_K0.1classified, samples_df$sample_type)

describeBy(samples_2$seq_K0.1classified, samples_2$collection_day)
ANOVA<-aov(samples_2$seq_K0.1classified~samples_2$collection_day)
summary(ANOVA)
kruskal.test(samples_2$seq_K0.1classified, samples_2$collection_day)
pairwise.t.test(samples_2$seq_K0.1classified,samples_2$collection_day, p.adj = "bonferroni")


```

# 3. Phyloseq object preparation
We will:
* Subset data: Samples (n=95), negative controls (n=3) and positive controls (n=2)
* Normalize (using the Cumulative Sum Scaling) the counts for all our samples to make the groups comparable between each other. We will normalize the negative and positive controls as well only to present that result.
* Our normalized samples will be subset to: (1) Bacteria and (2) Archaea. 

## 3.1. Subset raw data 
Samples (n=95), negative (n=3) and positive controls (n=2).
Samples Pre-weaning (n=32), At weaning (n=32) and Post-weaning (n=31)
Define your phyloseq to use (e.g. Rumen or Rumen0)
```{r}
# Change the right part of the equivalence 
  ps = Rumen
```

Obtain: Rumen_Mock, Rumen_Blank, Rumen_Sample, Rumen1, Rumen2, Rumen3
```{r}
  Rumen_Mock <- subset_samples(ps, sample_type =="Mock")
  Rumen_Blank <- subset_samples(ps,sample_type =="Blank") 
  Rumen_Sample <- subset_samples(ps, sample_type =="Sample")
  
  Rumen1 <- subset_samples(Rumen_Sample, collection_day =="Pre_weaning")
  Rumen2 <- subset_samples(Rumen_Sample, collection_day =="At_weaning")
  Rumen3 <- subset_samples(Rumen_Sample, collection_day =="Post_weaning")
```

## 3.2. Normalization (CSS)
Method Cumulative Sum Scaling (using metagenomeSeq package) 
Define the phyloseq object that you want to normalize: e.g. Rumen_Sample, Rumen_Blank, Rumen1, etc
```{r}
# Change the right part of the equivalence 
  ps = Rumen_Sample
```
NOTE: To obtain the number of classified reads at any given level, just agglomerate Rumen_Sample object to the level you want to investigate and then do: 
  sum(taxa_sums(glom)) # where glom is the agglomerated object

Obtain the CSS normalized object:
```{r}
# Changing from phyloseq object to MRexperiment (metagenomeSeq format) for CSS normalization
rumen.metaseq <- phyloseq_to_metagenomeSeq(ps) # change here

# CSS normalization (cumNorm)
# cumNormStat: Calculates the percentile for which to sum counts up to and scale by. 
rumen.metaseq.norm<- cumNorm(rumen.metaseq, p=cumNormStat(rumen.metaseq))

# Extracting (MRcounts) the normalized table to be used as otu_table (norm=TRUE)
CSS_rumen.metaseq <- MRcounts(rumen.metaseq.norm, norm = TRUE)

# Merging the extracted normalized count table to use as otu_table in normalized phyloseq object "CSS_Rumen"
CSS_Rumen <- merge_phyloseq(otu_table(CSS_rumen.metaseq,
taxa_are_rows=T),sample_data(ps),tax_table(ps))
```

## 3.3. Subset norm object
We will use *CSS_Rumen* normalized object obtained in the step before. We will subset different domains/kingdoms: 
- Bacteria: Rumen_bac
- Archaea: Rumen_arc 
- Virus: Rumen_vir
- Protozoa: Rumen_euk
- Fungi: Rumen_fun
```{r}
  Rumen_bac <- subset_taxa(CSS_Rumen, Domain %in% c("Bacteria"))
  Rumen_arc <- subset_taxa(CSS_Rumen, Domain %in% c("Archaea"))
  Rumen_vir <- subset_taxa(CSS_Rumen, Domain %in% c("Viruses"))
  Rumen_euk <- subset_taxa(CSS_Rumen, Domain %in% c("Eukaryota"))
  Rumen_pro <- subset_taxa(Rumen_euk, Kingdom %in% c("NA"))
  Rumen_fun <- subset_taxa(Rumen_euk, Kingdom %in% c("Fungi"))
```

For confidence parameter 0.1 (Rumen0)
*Rumen_bac:* 9819 unique OTUs (taxa)
*Rumen_arc:* 419 unique OTUs (taxa)
*Rumen_vir:* 495 unique OTUs (taxa)
*Rumen_pro:* 47 unique OTUs (taxa)
*Rumen_fun:* 112 unique OTUs (taxa)

## 3.4. Agglomerate taxa to different levels
Define the ps object that you will agglomerate (e.g., Rumen_bac, Rumen_arc, etc)
Also define the level you will agglomerate (e.g., Phylum, Genus, Species, etc)
```{r}
# Change the right part of the equivalence 
  ps = Rumen_fun
  lvl = "Species"
```
Agglomerate taxa:
```{r}
  glom= tax_glom(ps, lvl)
```
Recover object -- Give a proper name to the agglomerated object according to the kingdom/domain you agglomerate and the level you agglomerated (e.g., Rumen_bac_Phy, Rumen_arc_Phy, etc)
```{r}
# change the left part of the equivalence 
  Rumen_arc_Gen = glom
```

For confidence parameter 0.1 (Rumen0)
*Rumen_bac_Phy:* 40 unique Phylum
*Rumen_arc_Phy:* 7 unique Phylum
*Rumen_vir_Phy:* 14 unique Phylum
*Rumen_pro_Phy:* 6 unique Phylum
*Rumen_fun_Phy:* 3 unique Phylum

*Rumen_bac_Gen:* 1909 unique Genus
*Rumen_arc_Gen:* 140 unique Genus
*Rumen_vir_Gen:* 368 unique Genus
*Rumen_pro_Gen:* 18 unique Genus
*Rumen_fun_Gen:* 48 unique Genus

# 4. Relative abundance
We will analyze to Phylum and Genus level (and maybe other depending on the domain/kingdom)
## 4.1. Controls
```{r}
# Change the right part of the equivalence 
  ps = Rumen_Mock_Gen
```
Obtain the relative_long for ggplot:
```{r}
# 1. Converting to Percentage and transforming Genus tax_glom object to long-format table
  relative <- transform_sample_counts(ps, function(x) (x / sum(x))*100 )
  relative_long <- psmelt(relative)
# 2. Creating a new column (variable) for the  long-format table, called "mean_relative_abund" 
  relative_long <- relative_long %>%
  group_by(Genus) %>%
  mutate(mean_relative_abund = mean(Abundance))
# 3. Formatting the "Genus" column variable as names and the "mean_relative_abund" column variable as continuous numbers
  relative_long$Genus <- as.character(relative_long$Genus)
  relative_long$mean_relative_abund <- as.numeric(relative_long$mean_relative_abund)
# 4. Giving the name "Others (< 1%)" to all the OTUs in the "Genus" column that have less than 1 in the "Abundace" column in the long-format table 
  relative_long$Genus[relative_long$Abundance < 0.1] <- "Others (< 0.1%)" 
```
### SUP FIGURE 2. Controls
```{r}
# POSITIVE controls
SupF1 <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
# NEGATIVE controls
SupF2 <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
```

## 4.2. Phylum level
```{r}
# Change the right part of the equivalence 
  ps = Rumen_vir_Phy
```
Obtain the relative_long for ggplot:
```{r}
# 1. Converting to Percentage and transforming Phylum tax_glom object to long-format table
  relative <- transform_sample_counts(ps, function(x) (x / sum(x))*100 )
  relative_long <- psmelt(relative)
# 2. Creating a new column (variable) for the  long-format table, called "mean_relative_abund" 
  relative_long <- relative_long %>%
  group_by(Phylum) %>%
  mutate(mean_relative_abund = mean(Abundance))
# 3. Formatting the "Phylum" column variable as names and the "mean_relative_abund" column variable as continuous numbers
  relative_long$Phylum <- as.character(relative_long$Phylum)
  relative_long$mean_relative_abund <- as.numeric(relative_long$mean_relative_abund)
# 4. Giving the name "Others (< 1%)" to all the OTUs in the "Genus" column that have less than 1 in the "Abundace" column in the long-format table 
  relative_long$Phylum[relative_long$Abundance < 1] <- "Others (< 1%)" 
```
Exploratory figure: Relative abundance per collection day
```{r}
relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
```

### FIGURE 1. Classification at phylum level
```{r}
# A. Archaea
relative_long$Phylum = factor(relative_long$Phylum, levels=c("Others (< 1%)","Candidatus Thermoplasmatota","Euryarchaeota" )) # Change order in plot HERE
F1A <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=15),legend.title=element_text(size=17))+
  ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20") 


# B. Bacteria
relative_long$Phylum = factor(relative_long$Phylum, levels=c("Others (< 1%)","Spirochaetes", "Fibrobacteres","Actinomycetota","Pseudomonadota","Bacillota","Bacteroidota")) # Change order in plot HERE

F1B <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=15),legend.title=element_text(size=17))+
  ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")

# C. Fungi
relative_long$Phylum = factor(relative_long$Phylum, levels=c("Others (< 1%)","Microsporidia","Basidiomycota","Ascomycota")) # Change order in plot HERE
F1C <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=15),legend.title=element_text(size=17))+
  ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")

# D. Protozoa
relative_long$Phylum = factor(relative_long$Phylum, levels=c("Others (< 1%)","Heterolobosea","Fornicata","Bacillariophyta", "Evosea",  "Euglenozoa", "Apicomplexa")) # Change order in plot HERE
F1D <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=15),legend.title=element_text(size=17))+
  ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")

# E. Virus
relative_long$Phylum = factor(relative_long$Phylum, levels=c("Others (< 1%)", "Artverviricota","Cossaviricota", "Cressdnaviricota", "Dividoviricota","Duplornaviricota","Kitrinoviricota","NA", "Negarnaviricota","Nucleocytoviricota","Peploviricota","Pisuviricota","Preplasmiviricota",     "Uroviricota")) # Change order in plot HERE
F1E <- relative_long %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank(),legend.text=element_text(size=15),legend.title=element_text(size=17))+
  ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
```
Obtain phylum-level relative abundance for Bacteria
```{r}
# Overall Phyla abundance in bacteria
  bac_Phy<-table1(~ Abundance | Phylum, data=relative_long)
```

## 4.3. Genus level
```{r}
# Change the right part of the equivalence 
  ps = Rumen_bac_Gen
```
Obtain the relative_long for ggplot:
```{r}
# 1. Converting to Percentage and transforming Genus tax_glom object to long-format table
  relative2 <- transform_sample_counts(ps, function(x) (x / sum(x))*100 )
  relative_long2 <- psmelt(relative2)
# 2. Creating a new column (variable) for the  long-format table, called "mean_relative_abund" 
  relative_long2 <- relative_long2 %>%
  group_by(Genus) %>%
  mutate(mean_relative_abund = mean(Abundance))
# 3. Formatting the "Genus" column variable as names and the "mean_relative_abund" column variable as continuous numbers
  relative_long2$Genus <- as.character(relative_long2$Genus)
  relative_long2$mean_relative_abund <- as.numeric(relative_long2$mean_relative_abund)
# 4. Giving the name "Others (< 1%)" to all the OTUs in the "Genus" column that have less than 1 in the "Abundace" column in the long-format table 
  relative_long2$Genus[relative_long2$Abundance < 1] <- "Others (< 1%)" 
```

###SUP FIGURE 3. Classification at genus level
```{r}
# BACTERIA
SupF3A <- relative_long2 %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")

# ARCHAEA
SupF3B <- relative_long2 %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_bar(stat = "identity", alpha=0.9)+theme_bw()+
  facet_wrap(~collection_day+weaning, nrow=1, scales="free_x")+
  theme(axis.text.x = element_blank())+ylab("Relative abundance (%)")+
  ggthemes::scale_fill_tableau("Tableau 20")
```
Obtain tables for relative abundance
```{r}
 # Overall Genera abundance in bacteria
  arc_Gen <- table1(~ Abundance | Genus, data=relative_long2)
  bac_Gen <- table1(~ Abundance | Genus, data=relative_long2)
```

# 5. Alpha diversity
We will estimate the alpha diversity indices and plot them. Finally we will modelate properly our indices as outcomes and the predictors of our study.
## 5.1. Estimate alpha diversity
Define your raw phyloseq object to be subsetted and the Domains you want to subset. Also define your level for taxa_glom (Usually genus)
```{r}
# Change the right part of the equivalence 
  ps = Rumen
  div = c("Bacteria", "Archaea") # Domain(s) required
  lvl = "Genus"
```
Subset the domain(s) and agglomerate for alpha diversity
```{r}
  ps <- subset_taxa(ps, Domain %in% c(div))
  glom= tax_glom(ps, lvl)
```
Recover object -- Give a proper name to the agglomerated object according to the the level you agglomerated
```{r}
# change the left part of the equivalence 
  Rumen_div = glom
```
Estimate the richness, Simpson, Shannon and Evenness indices
```{r}
  div1= estimate_richness(Rumen_div, measures=c("Observed", "InvSimpson", "Shannon"))
  div2= evenness(Rumen_div, 'pielou')
```
Save the results in a table
```{r}
# Change the name of the csv file. div1 is for richness, simpson, shannon) and div2 is for evenness
  write.csv(div1, "Rumen_div1_final.csv")
  write.csv(div2, "Rumen_div2_final.csv")
```

## 5.2. Plotting alpha diversity
Re-level the categorical variables
```{r}
samples_2$castration_group=factor(samples_2$castration_group, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning"))
samples_2$castration=factor(samples_2$castration, levels=c("Birth", "Turnout", "Pre_weaning", "Weaning", "Not_castrated"))
samples_2$weaning_group=factor(samples_2$weaning_group, levels=c("Fence_line", "Truck"))
samples_2$weaning=factor(samples_2$weaning, levels=c("Fence_line", "Truck", "Not_weaned"))
samples_2$collection_day=factor(samples_2$collection_day, levels=c("Pre_weaning", "At_weaning", "Post_weaning"))
```

Define the variables for the plot:
y_axis can be   = Kra0.1_Obs_Genus, Kra0.1_Shan_Genus, Kra0.1_Simp_Genus, Kra0.1_Even_Genus
y_label can be  = Richness / Shannon Index / Inv-Simpson Index / Pielou's evenness index

FOR colorby = collection_day 
  color = c("brown", "#D98326", "#E8B600")
FOR colorby = weaning_group or weaning
  color = c("#967bba", "#505359", "#c9c8c8")
FOR colorby = castration_group / castration
  color = c("#f54545", "#85a9ca", "#e0b54e", "#73a054", "#c9c8c8")
```{r}
# Change the right part of the formula
  x_axis = "collection_day"
  x_label = "Collection day"
  y_axis = "Kra0.1_Shan_Genus"
  y_label = "Shannon Index"
  colorby = "weaning_group"
  color = c("#967bba", "#505359", "#c9c8c8")
```
Exploratory figure. Note that you will need to change the dataframe "samples_2" according to your needs
```{r}
# plotting richness, shannon and Pielou by WEANING and Collection day
  ggplot(samples_2,
         aes_string(x_axis,y_axis, color=colorby))+
  geom_boxplot(position=position_dodge(0.8),lwd=1)+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+
  theme_classic()+ 
  #ylim(0,1)+
  geom_jitter(position=position_dodge(0.8),size=2)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  #theme(axis.text.x = element_blank())+
  ylab(y_label)+
  xlab(x_label)
```
## FIGURE 2B
```{r}
F2B1 <-   ggplot(samples_2,
         aes_string(x_axis,y_axis, color=colorby))+
  geom_boxplot(position=position_dodge(0.8),lwd=1)+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+
  theme_classic()+ 
  #ylim(0,1)+
  geom_jitter(position=position_dodge(0.8),size=2)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  ylab(y_label)+
  xlab(x_label)
```
## FIGURE 3B
```{r}
F3B <-   ggplot(samples_2,
         aes_string(x_axis,y_axis, color=colorby))+
  geom_boxplot(position=position_dodge(0.8),lwd=1)+
  scale_fill_manual(values = color)+
  scale_color_manual(values = color)+
  theme_classic()+ 
  #ylim(0,1)+
  geom_jitter(position=position_dodge(0.8),size=2)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  ylab(y_label)+
  xlab(x_label)
```

## 5.3. Statistics 
We used a linear mixed-effect model to fit random (cow_ID) and fixed effects (castration / castration_group, weaning / weaning_group, age_days, weight_ADGUpToDate)
* Our outcome (response variable) will be genus-level alpha indices: richness and shannon 
* Our predictor (explanatory variable) will be the interventions and time: castration, weaning and collection_days
```{r}
# Define data frames
md= samples_2
md1= samples_D1
md2= samples_D2
md3= samples_D3
```
Run model for WEANING
```{r}
  weaning_shan = lmer(Kra0.1_Shan_Genus ~ weaning_group*collection_day +castration_group + age_days + weight_ADGUpToDate+ (1|cow_ID), data=md, REML = T)
  # weaning_shan1 = lmer(Kra0.1_Shan_Genus ~ weaning_group*collection_day +castration_group + age_days + weight_ADGUpToDate+ (1|cow_ID), data=md, REML = F)
  # weaning_shan2 = lmer(Kra0.1_Shan_Genus ~ weaning_group*collection_day +castration + age_days + weight_ADGUpToDate+ (1|cow_ID), data=md, REML = F)  
  # anova(weaning_shan1, weaning_shan2 )
  summary(weaning_shan)
  confint(weaning_shan)
  emmeans(weaning_shan,pairwise~weaning_group|collection_day)
  Anova(weaning_shan, type = "III")
  Anova(weaning_shan, type = "III", test.statistic = "F")

  weaning_obs = lmer(Kra0.1_Obs_Genus ~ weaning_group*collection_day +castration_group + age_days + weight_ADGUpToDate+ (1|cow_ID), data=md, REML = T)
  summary(weaning_obs)
  confint(weaning_obs)
  emmeans(weaning_obs,pairwise~weaning_group|collection_day)
  Anova(weaning_obs, type = "III")
  Anova(weaning_obs, type = "III", test.statistic = "F")
```
Run model for CASTRATION
```{r}
# CASTRATION
  castration_shan = lmer(Kra0.1_Shan_Genus~castration_group*collection_day+weaning_group + (1|cow_ID), data=md, REML = T)
  # anova(castration_shan1, castration_shan2)
  summary(castration_shan)
  emmeans(castration_shan,pairwise~castration_group|collection_day)
  Anova(castration_shan, type = "III")
  Anova(castration_shan, type = "III", test.statistic = "F")
  
  castration_obs = lmer(Kra0.1_Obs_Genus~castration_group*collection_day+weaning_group + (1|cow_ID), data=md, REML = T)
  summary(castration_obs)
  emmeans(castration_obs,pairwise~castration_group|collection_day)
  Anova(castration_obs, type = "III")
  Anova(castration_obs, type = "III", test.statistic = "F")
  
# CASTRATION PER DAYS
  castration_shan1 = lm(Kra0.1_Shan_Genus~castration + age_days + weight_ADGUpToDate, data=md1)
  summary(castration_shan1)
  emmeans(castration_shan1,pairwise~castration) 
  Anova(castration_shan1, type = "III")
 
  castration_shan2 = lm(Kra0.1_Shan_Genus~castration + age_days + weight_ADGUpToDate, data=md2)
  summary(castration_shan2)
  emmeans(castration_shan2,pairwise~castration) 
  Anova(castration_shan2, type = "III")
  
  castration_shan3 = lm(Kra0.1_Shan_Genus~castration*weaning_group + age_days + weight_ADGUpToDate, data=md3)
  summary(castration_shan3)
  emmeans(castration_shan3,pairwise~castration) 
  Anova(castration_shan3, type = "III")
  
# WEANING PER DAY
    weaning_shan1 = lm(Kra0.1_Shan_Genus~weaning_group + age_days + weight_ADGUpToDate, data=md1)
    weaning_shan2 = lm(Kra0.1_Shan_Genus~weaning_group + age_days + weight_ADGUpToDate, data=md2)
    weaning_shan3 = lm(Kra0.1_Shan_Genus~weaning_group + age_days + weight_ADGUpToDate, data=md3)
    summary(weaning_shan3)
```

# 6. Beta diversity
We will estimate the *ordination* for the plot and the *distance matrix* for the statistics
## 6.1. Estimate beta diversity
Define your CSS-normalized phyloseq object to be subsetted (e.g., CSS_Rumen, CSS_Rumen1, CSS_Rumen2, CSS_Rumen3) and the Domains you want to subset. Also define your level for taxa_glom (Usually genus)
```{r}
# Change the right part of the equivalence 
  ps = CSS_Rumen
  div = c("Bacteria", "Archaea") # Domain(s) required
  lvl = "Genus"
```
Subset the domain(s) and agglomerate for alpha diversity
```{r}
  ps <- subset_taxa(ps, Domain %in% c(div))
  glom = tax_glom(ps, lvl)
```
Recover object -- Give a proper name to the agglomerated object according to the the level you agglomerated (e.g., Rumen_div_Gen, Rumen_div_Gen1, etc)
```{r}
# change the left part of the equivalence 
  Rumen_div_Gen = glom
```

Define your CSS-normalized genus-agglomerate phyloseq object from which we will get the ordination.
```{r}
# Change the right part of the equivalence 
  ps = Rumen_div_Gen_clean
  set.seed(143) # it will set a fixed number for all the random-based permutational tests, so it can be reproducible
```
Get the ordination using NMDS and bray curtis distance matrix. FOR PLOT
```{r}
  ord <- ordinate(ps, "NMDS", "bray")
```
Recover object -- Give a proper name to the ordination (e.g., Rumen_ord, Rumen_ord1, etc)
```{r}
# change the left part of the equivalence 
  Rumen_ord_clean = ord
```

Define your CSS-normalized genus-agglomerate phyloseq object from which we will get the distance matrix (dist)
```{r}
# Change the right part of the equivalence 
  ps = Rumen_div_Gen_clean
# Rumen_div_Gen_clean = removed R0048 because did not allow to get ordination for collection_day 1 and improved the figure but DID NOT CHANGE THE DIST MATRIX NOR STATISTICS 
```
Get the distance matrix. FOR STATISTICS (permanova)
```{r}
  dist = vegdist(t(otu_table(ps)), method ="bray" )
```
Recover object -- Give a proper name to the distance matrix (e.g., Rumen_dist, Rumen_dist1, etc)
```{r}
# change the left part of the equivalence 
  Rumen_dist_clean = dist
# Rumen_dist_clean = removed R0048 because did not allow to get ordination for collection_day 1 and improved the figure but DID NOT CHANGE THE DIST MATRIX NOR STATISTICS 
```

## 6.2. Plotting beta diversity
Define the variables for the plot:
FOR colorby = collection_day 
  color = c("brown", "#D98326", "#E8B600")
FOR colorby = weaning_group or weaning
  color = c("#967bba", "#505359", "#c9c8c8")
FOR colorby = castration_group / castration
  color = c("#f54545", "#85a9ca", "#e0b54e", "#73a054", "#c9c8c8")
FOR shapeby = collection_day
  shape = c(15,16,17)
```{r}
# Change the right part of the formula
  ps = Rumen_div_Gen_clean
  ord = Rumen_ord_clean
  colorby = "weaning"
  color = c("#967bba", "#505359", "#c9c8c8")
  shapeby = "collection_day" # OPTIONAL
  shape = c(15,16,17) # OPTIONAL
```
Exploratory figure. Note that you will need to change the dataframe "samples_2" according to your needs
```{r}
p<-plot_ordination(ps, ord, type="sample", color=colorby,
                   #shape="collection_day")+
                   )+
  geom_point(size=4.5, alpha=0.8)+
  scale_fill_manual(values = color)+
  #scale_shape_manual(values = shape)+
  scale_color_manual(values = color)+
  stat_ellipse()+
  theme_classic() 

p$layers
p$layers <- p$layers[-1]
p
```
## FIGURE 2A
```{r}
F2A <-plot_ordination(ps, ord, type="sample", color=colorby,
                   #shape="collection_day")+
                   )+
  geom_point(size=4.5, alpha=0.8)+
  scale_fill_manual(values = color)+
  #scale_shape_manual(values = shape)+
  scale_color_manual(values = color)+
  stat_ellipse()+
  theme_classic() 

F2A$layers
F2A$layers <- F2A$layers[-1]
F2A
```
## FIGURE 3A
```{r}
levels(sample_data(ps)$weaning) <- c("Fence_line", "Truck","Before_weaning") # Re-labeling data

F3A <-plot_ordination(ps, ord, type="sample", color=colorby,
                   shape="collection_day")+
  geom_point(size=4.5, alpha=0.8)+
  scale_fill_manual(values = color)+
  scale_shape_manual(values = shape)+
  scale_color_manual(values = color)+
  stat_ellipse()+
  theme_classic() 

F3A$layers
F3A$layers <- F3A$layers[-1]
F3A
```
## FIGURE 4B
```{r}
F4B1.2 <-plot_ordination(ps, ord, type="sample", color=colorby,
                   #shape="collection_day")+
                   )+
  geom_point(size=4.5, alpha=0.8)+
  scale_fill_manual(values = color)+
  #scale_shape_manual(values = shape)+
  scale_color_manual(values = color)+
  stat_ellipse()+
  theme_classic() 

F4B1.2$layers
F4B1.2$layers <- F4B1.2$layers[-1]
F4B1.2
```

## 6.3. Statistics
We used a use PERMANOVA (adonis2) to determine the variation explained by our predictors (castration / castration_group, weaning / weaning_group, collection_day)
* Our outcome (response variable) will be distance matrix 
* Our predictor (explanatory variable) will be the interventions and time: castration, weaning and collection_days

Define data frames (metadata)
```{r}
md # All
md1 # collection_day1
md2 # collection_day2
md3 # collection_day3
md_clean= data.frame(sample_data(Rumen_div_Gen_clean))
# md_clean = removed R0048 because did not allow to get ordination for collection_day 1 and improved the figure but DID NOT CHANGE THE DIST MATRIX NOR STATISTICS 
md1_clean= data.frame(sample_data(Rumen_div_Gen1.2))

set.seed(143) # it will set a fixed number for all the random-based permutational tests, so it can be reproducible
```
Run model for WEANING
```{r}
weaning_beta = adonis2(Rumen_dist_clean ~ weaning_group  + collection_day + castration, data=md_clean , by="margin")
# Rumen_dist_clean = removed R0048 because did not allow to get ordination for collection_day 1 and improved the figure but DID NOT CHANGE THE DIST MATRIX NOR STATISTICS 

# Per day
  weaning_beta1 = adonis2(Rumen_dist1.2 ~ weaning_group + castration, data=md1_clean)
  weaning_beta2 = adonis2(Rumen_dist2 ~ weaning_group + castration, data=md2)
  weaning_beta3 = adonis2(Rumen_dist3 ~ weaning_group + castration, data=md3)
```
Run model for CASTRATION for each collection_day
```{r}
  castration_beta1 = adonis2(Rumen_dist1.2 ~ castration , data=md1_clean) 
  castration_beta2 = adonis2(Rumen_dist2 ~ castration , data=md2) 
  castration_beta3 = adonis2(Rumen_dist3 ~ castration + weaning_group , data=md3) # added weaning as confounder because the castration groups were weaned as well
```
## Extra: Beta-dispersion
```{r}
## BETA DISPERSION

## DISPERSION: Check this later ##

disp.weaning = betadisper(dist, type = c("centroid"), metadata$collection_day)
anova(disp.weaning) # ns
TukeyHSD(disp.weaning, which = "group", ordered = FALSE,
         conf.level = 0.95)

#permutest(disp.weaning, pairwise=TRUE, permutations=999)
boxplot(disp.weaning)
```

# 7. Differential abundance
We are going to account for non-normal distribution of OTU's and sparse (many 0) of microbiome data using a "Zero-inflated Gaussian Mixture Model" to assess statistical significance on microbial differential abundance between collection days (time) and, castration and weaning strategies. Pairwise comparisons will be done and we will use the adjusted p-value (BH correction) for multiple comparisons and Log2Fold Change to decide which microbes where significant both statistically (p-value) and biologically (logFC) different on abundance over groups of contrast. 
We will focus for now on:
* collection_day at phylum level (Contrasting: At_weaning vs Pre_weaning; Post_weaning vs Pre_weaning)
* weaning_group at genus level (Contrasting: Fence_line vs Truck IN Pre_weaning; Fence_line vs Truck IN At_weaning; Fence_line vs Truck IN Post_weaning)
* castration at phylum level (Constrasting: 
IN Pre_weaning: Birth vs Not_castrated, Turnout vs Not_castrated, Turnout vs Birth;
IN At_weaning: Birth vs Not_castrated, Turnout vs Not_castrated, Pre_weaning vs Not_castrated, Pre_weaning vs Turnout, Pre_weaning vs Birth, Turnout vs Birth;
IN Post_weaning: Weaning vs Pre_weaning,  Weaning vs Turnout, Weaning vs Birth, Pre_weaning vs Turnout, Pre_weaning vs Birth; Turnout vs Birth)

The strings for all this comparisons are here:
for contrasting you can have any of the following:
c("collection_dayAt_weaning-collection_dayPre_weaning","collection_dayPost_weaning-collection_dayPre_weaning")
c("weaning_groupFence_line-weaning_groupTruck")
c("castrationBirth-castrationNot_castrated","castrationTurnout-castrationNot_castrated","castrationTurnout-castrationBirth")
c("castrationBirth-castrationNot_castrated","castrationTurnout-castrationNot_castrated", "castrationPre_weaning-castrationNot_castrated","castrationPre_weaning-castrationTurnout","castrationPre_weaning-castrationBirth","castrationTurnout-castrationBirth")
c("castrationWeaning-castrationPre_weaning","castrationWeaning-castrationTurnout","castrationWeaning-castrationBirth","castrationPre_weaning-castrationTurnout","castrationPre_weaning-castrationBirth","castrationTurnout-castrationBirth")

## 7.1. collection_day at phylum level
Define your variables for ACROSS TIME CODE. This will include random effects of repetitive measures over the same animal.
```{r}
# change the right part of the formula
  ps = CSS_Rumen # object to be analyzed (it can be: CSS_Rumen, CSS_Rumen1, CSS_Rumen2, CSS_Rumen3)
  lvl = "Phylum" # level of the differential abundance testing (It can be: Phylum, Genus)
  raw = Rumen_Sample # for model zero that accounts for depth of coverage (It can be Rumen, Rumen1, Rumen2, Rumen3)
  contrasting = c("collection_dayAt_weaning-collection_dayPre_weaning","collection_dayPost_weaning-collection_dayPre_weaning") # depends on what comparison are you interested
```
Run the *ACROSS TIME CODE (random effects) - Part 1*
Code for correlated data, meaning the observation are non-independent (Independence assumption is violated)
```{r}
# 1. DEFINE YOUR METAGENOMESEQ EXPERIMENT (cumNorm-ed)
# we need to filter the OTUs at Phylum level to be present in at least in 2 samples (see ?filterData for help)
  metaseq <- phyloseq_to_metagenomeSeq(ps) # it will change to MRexperiment as required by the "filterData" function
  filtered <- filterData(aggTax(metaseq, lvl=lvl), present = 2) # it will aggregate taxa at phylum (and give phylum names instead of only "OTU XXXX"), and filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  filtered.metaseq <- cumNorm(filtered) # it will normalize again the new filtered MRexperiment 
# 2. DEFINE YOUR METADATA VARIABLES TO INCLUDE IN THE MODEL (pData)
  collection_day <- pData(filtered.metaseq)$collection_day
  castration <-  pData(filtered.metaseq)$castration # will include 4 castration strategies + Not_castrated
  weaning_group <-  pData(filtered.metaseq)$weaning_group # will include 2 weaning strategies 
  cow_ID <- pData(filtered.metaseq)$cow_ID
# 3. Define your settings (zigControl), we will use the ones by default
  settings = zigControl(maxit = 10, verbose = TRUE)  
# 4. DEFINE YOUR MODEL ZERO TO ACCOUNT FOR DEPTH OF COVERAGE (zero_mod)
## MODEL 0: ZERO MODEL FOR ALL DAYS (depth of coverage of the NOT-normalized (RAW) MRexperiment)
  raw.metaseq <- phyloseq_to_metagenomeSeq(raw) # change here
  raw.filtered.metaseq <- filterData(aggTax(raw.metaseq,lvl=lvl), present = 2) # it will filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  zero_mod <- model.matrix(~0+log(libSize(raw.filtered.metaseq))) # model it!
# 5. DEFINE YOUR MODEL TO TEST (model.matrix), FIT THE MODEL (fitZig), MAKE PAIRWISE CONTRASTS (makeContrasts, contrasts.fit, eBayes) & EXTRACT RESULTS (topTable) 
# MODEL 1: Collection day
# Assign model (model.matrix)
  design_collection_day = model.matrix(~0 + collection_day)
  design_weaning = model.matrix(~0 +weaning_group)
  design_castration = model.matrix(~0 +castration) 
```

Define your model
```{r}
  model = design_collection_day # define the variable of contrast (it can be: design_weaning, design_castration, design_collection_day)
```
Run the *ACROSS TIME CODE (random effects) - Part 2*
```{r}
# Assign correlated data (random effect)
  voom_model <- voom(MRcounts(filtered.metaseq), model)
  dup <- duplicateCorrelation(voom_model, block=cow_ID, design= model)
  dup$consensus.correlation # to check your correlation coeff?
## Check Noelle's question for duplicateCorrelation syntax:
  #https://support.bioconductor.org/p/77336/
# Make ZIG model 
  ZIG = fitZig(obj= filtered.metaseq, mod = model, control = settings,zeroMod=zero_mod, useCSSoffset=FALSE, useMixedModel=dup$consensus.correlation, block=cow_ID)
  ZIGfit = ZIG@fit
  Mod = ZIG@fit$design
# Contrasts by collection dayp
  contrast = makeContrasts(contrasts=contrasting, levels=Mod)
  res = contrasts.fit(ZIGfit, contrast)
  resEB = eBayes(res)
```

Define the variables for the csv file with the table you want to recover
```{r}
# change the right part of the formula
  numb = 1 # number of comparison you want to recover according to your "contrasting" string
  tablename = "DALL_collectionday_AW-PreW.csv" # name the csv file according to the comparison you did (e.g. "D1_weaning_FL-T.csv", "D2_castration_TO-B.csv")  
```
Recover the table
```{r}
# FINALLY EXTRACT YOUR RESULTS
  fz_table <- topTable(resEB, adjust.method="BH", coef = numb,number = 1000000)
  write.csv(fz_table, tablename)
```

## 7.2. weaning_group at genus level
Define your variables for SINGLE TIME-POINT CODE. This will NOT include random effects because you have only one measure of the animals.
```{r}
# change the right part of the formula
  ps = CSS_Rumen1 # object to be analyzed (it can be: CSS_Rumen, CSS_Rumen1, CSS_Rumen2, CSS_Rumen3)
  lvl = "Genus" # level of the differential abundance testing (It can be: Phylum, Genus)
  raw = Rumen1 # for model zero that accounts for depth of coverage (It can be Rumen, Rumen1, Rumen2, Rumen3)
  day = "Pre_weaning" # day you want to analyze (It can be: Pre_weaning, At_weaning, Post_weaning)
  contrasting = c("weaning_groupFence_line-weaning_groupTruck") # depends on what comparison are you interested
```
Run the *SINGLE TIME-POINT CODE (no random effects)*
Code for non-correlated data, meaning the observation are independent (Independence assumption of data holds)
```{r}
# 1. DEFINE YOUR METAGENOMESEQ EXPERIMENT (cumNorm-ed)
# we need to filter the OTUs at Phylum level to be present in at least in 2 samples (see ?filterData for help)
  metaseq <- phyloseq_to_metagenomeSeq(ps) # it will change to MRexperiment as required by the "filterData" function
  filtered <- filterData(aggTax(metaseq, lvl=lvl), present = 2) # it will aggregate taxa at phylum (and give phylum names instead of only "OTU XXXX"), and filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  filtered.metaseq <- cumNorm(filtered) # it will normalize again the new filtered MRexperiment 
# 2. DEFINE YOUR METADATA VARIABLES TO INCLUDE IN THE MODEL (pData)
  weaning_group <-  pData(filtered.metaseq)$weaning_group
  castration <-  pData(filtered.metaseq)$castration
# 3. Define your settings (zigControl), we will use the ones by default
  settings = zigControl(maxit = 10, verbose = TRUE)
# 4. DEFINE YOUR MODEL ZERO TO ACCOUNT FOR DEPTH OF COVERAGE (zero_mod)
## MODEL 0: ZERO MODEL FOR ALL DAYS (depth of coverage of the NOT-normalized (RAW) MRexperiment)
  raw.metaseq <- phyloseq_to_metagenomeSeq(raw) # change here
  raw.filtered.metaseq <- filterData(aggTax(raw.metaseq,lvl=lvl), present = 2) # it will filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  zero_mod <- model.matrix(~0+log(libSize(raw.filtered.metaseq))) # model it!
## MODEL 0: ZERO MODEL FOR DAY3
### Subset
  samplesToKeep=which(pData(raw.filtered.metaseq)$collection_day==day) 
  raw.filtered.metaseq_filtered=raw.filtered.metaseq[,samplesToKeep] 
  zero_mod <- model.matrix(~0+log(libSize(raw.filtered.metaseq_filtered)))
# 5. DEFINE YOUR MODEL TO TEST (model.matrix), FIT THE MODEL (fitZig), MAKE PAIRWISE CONTRASTS (makeContrasts, contrasts.fit, eBayes) & EXTRACT RESULTS (topTable) 
# MODEL 1: Castration DAY 3
# Assign model (model.matrix)
  design_weaning = model.matrix(~0 + weaning_group)
  design_castration = model.matrix(~0 + castration)
```

Define your model
```{r}
  model = design_weaning # define the variable of contrast (it can be: design_weaning, design_castration, design_collection_day)
```
Run the *SINGLE TIME-POINT CODE (no random effects) - Part 1*
```{r}
# Make ZIG model 
  ZIG = fitZig(obj= filtered.metaseq, mod = model, control = settings, zeroMod=zero_mod, useCSSoffset=FALSE)
  ZIGfit = ZIG@fit
  Mod = ZIG@fit$design
# Contrasts by castration group
  contrast= makeContrasts(contrasts=contrasting,levels=Mod)
  res = contrasts.fit(ZIGfit, contrast)
  resEB = eBayes(res)
```

Define the variables for the csv file with the table you want to recover
```{r}
# change the right part of the formula
  numb = 1 # number of comparison you want to recover according to your "contrasting" string
  tablename = "D1_weaning_FL-T.csv" # name the csv file according to the comparison you did (e.g. "D1_weaning_FL-T.csv", "D1_weaning_FL-T.csv")  
```
Recover the table
```{r}
# FINALLY EXTRACT YOUR RESULTS
  fz_table <- topTable(resEB, adjust.method="BH", coef = numb,number = 1000000)
  write.csv(fz_table, tablename)
```

## 7.3. castration at phylum level 
Define your variables for SINGLE TIME-POINT CODE. This will NOT include random effects because you have only one measure of the animals.
```{r}
# change the right part of the formula
  ps = CSS_Rumen1 # object to be analyzed (it can be: CSS_Rumen, CSS_Rumen1, CSS_Rumen2, CSS_Rumen3)
  lvl = "Phylum" # level of the differential abundance testing (It can be: Phylum, Genus)
  raw = Rumen1 # for model zero that accounts for depth of coverage (It can be Rumen, Rumen1, Rumen2, Rumen3)
  day = "Pre_weaning" # day you want to analyze (It can be: Pre_weaning, At_weaning, Post_weaning)
  contrasting = c("castrationBirth-castrationNot_castrated","castrationTurnout-castrationNot_castrated","castrationTurnout-castrationBirth")
```

Run the *SINGLE TIME-POINT CODE (no random effects) - Part 1*
Code for non-correlated data, meaning the observation are independent (Independence assumption of data holds)
```{r}
# 1. DEFINE YOUR METAGENOMESEQ EXPERIMENT (cumNorm-ed)
# we need to filter the OTUs at Phylum level to be present in at least in 2 samples (see ?filterData for help)
  metaseq <- phyloseq_to_metagenomeSeq(ps) # it will change to MRexperiment as required by the "filterData" function
  filtered <- filterData(aggTax(metaseq, lvl=lvl), present = 2) # it will aggregate taxa at phylum (and give phylum names instead of only "OTU XXXX"), and filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  filtered.metaseq <- cumNorm(filtered) # it will normalize again the new filtered MRexperiment 
# 2. DEFINE YOUR METADATA VARIABLES TO INCLUDE IN THE MODEL (pData)
  weaning_group <-  pData(filtered.metaseq)$weaning_group
  castration <-  pData(filtered.metaseq)$castration
# 3. Define your settings (zigControl), we will use the ones by default
  settings = zigControl(maxit = 10, verbose = TRUE)
# 4. DEFINE YOUR MODEL ZERO TO ACCOUNT FOR DEPTH OF COVERAGE (zero_mod)
## MODEL 0: ZERO MODEL FOR ALL DAYS (depth of coverage of the NOT-normalized (RAW) MRexperiment)
  raw.metaseq <- phyloseq_to_metagenomeSeq(raw) # change here
  raw.filtered.metaseq <- filterData(aggTax(raw.metaseq,lvl=lvl), present = 2) # it will filter the data from the MRexperiment. YOU CAN FILTER AS YOU WANT, TRY FILTERING BY DEPTH: "filterData(obj, present = X, depth = y)
  zero_mod <- model.matrix(~0+log(libSize(raw.filtered.metaseq))) # model it!
## MODEL 0: ZERO MODEL FOR DAY3
### Subset
  samplesToKeep=which(pData(raw.filtered.metaseq)$collection_day==day) 
  raw.filtered.metaseq_filtered=raw.filtered.metaseq[,samplesToKeep] 
  zero_mod <- model.matrix(~0+log(libSize(raw.filtered.metaseq_filtered)))
# 5. DEFINE YOUR MODEL TO TEST (model.matrix), FIT THE MODEL (fitZig), MAKE PAIRWISE CONTRASTS (makeContrasts, contrasts.fit, eBayes) & EXTRACT RESULTS (topTable) 
# MODEL 1: Castration DAY 3
# Assign model (model.matrix)
  design_weaning = model.matrix(~0 + weaning_group)
  design_castration = model.matrix(~0 + castration)
```

Define your model
```{r}
  model = design_castration # define the variable of contrast (it can be: design_weaning, design_castration, design_collection_day)
```
Run the *SINGLE TIME-POINT CODE (no random effects) - Part 2*
```{r}
# Make ZIG model 
  ZIG = fitZig(obj= filtered.metaseq, mod = model, control = settings, zeroMod=zero_mod, useCSSoffset=FALSE)
  ZIGfit = ZIG@fit
  Mod = ZIG@fit$design
# Contrasts by castration group
  contrast= makeContrasts(contrasts=contrasting,levels=Mod)
  res = contrasts.fit(ZIGfit, contrast)
  resEB = eBayes(res)
```

Define the variables for the csv file with the table you want to recover
```{r}
# change the right part of the formula
  numb = 1 # number of comparison you want to recover according to your "contrasting" string
  tablename = "D1_castration_B-NC.csv" # name the csv file according to the comparison you did (e.g. "D1_weaning_FL-T.csv", "D2_castration_TO-B.csv")  
```
Recover the table
```{r}
# FINALLY EXTRACT YOUR RESULTS
  fz_table <- topTable(resEB, adjust.method="BH", coef = numb,number = 1000000)
  write.csv(fz_table, tablename)
```

## 7.3. Volcano Plots

Small number of features: Phylum-level plots

Define you variables
```{r}
  Log_Phy <- read.csv("LogFC_Rumen_Phy.csv")
# Change the right hand part of the formula
  tittle = "Collection days: Phylum-level pairwise comparison"

  da = "All"
  vari = "collection_day"
  relevel = c("At weaning vs Pre-weaning", "Post-weaning vs At weaning")

# relevel can be:
# collection_day = c("At weaning vs Pre-weaning", "Post-weaning vs At weaning")
# PW_castration = c("Birth vs Not castrated", "Turn-out vs Birth", "Turn-out vs Not castrated")
# AW_castration = c("Birth vs Not castrated", "Pre-weaning vs Birth", "Pre-weaning vs Not castrated", "Pre-weaning vs Turn-out", "Turn-out vs Birth", "Turn-out vs Not castrated")
# PW_castration = c("Pre-weaning vs Birth", "Pre-weaning vs Turn-out", "Turn-out vs Birth", "Weaning vs Birth", "Weaning vs Pre-weaning", "Weaning vs Turn-out")
```
get your figure
```{r}
  # Subsetting by categorical predictor (variable) and getting rid of NAs and not_match
  tmp <- subset(Log_Phy, Log_Phy$variable == vari)
  tmp <- subset(tmp, tmp$day == da)
  tmp <- filter(tmp, taxa != "NA", taxa !="no_match")
#median(log2.phy$AveExpr)
  # re-leveling the levels of categorical variable
  tmp$comparison_type = factor(tmp$comparison_type, levels=relevel)
  # creating a new column (variable) for significant p-value
  tmp$Significant <- ifelse(tmp$adj.P.Val < 0.05, "AdjP < 0.05", "Not Significant")
  
# OTHER WAY TO FILTER DATA
#Log_Phy$EffectSize <- ifelse(Log_Phy$AveExpr < 4.2057385, "ES<1", "EF>1")
#topTable_contrasts_df['taxa'] <- row.names(topTable_contrasts_df)
#Log_Phy_fil <- subset(Log_Phy,  percentile_fifty_cut%in% c("1"))       ## taxa av.exp >50th percentile
#Log_Phy_fil
#size_AveEx <- sqrt(Log_Phy_fil$AveExpr/pi)    ## size of average expression/effectsize

  name <- ggplot(tmp, aes(x = logFC, y = taxa, fill = Significant))+
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
  ggtitle(tittle)+         # change the title based on contrast 1,2..
  theme(legend.position = 'none')
```

Recover the name of figure
```{r}
# change the left part of the formula
  F2C <- name 
```

Large number of features: Genus-level plots

Define you variables
```{r}
  Log_Gen <- read.csv("LogFC_Rumen_Gen.csv")
# Change the right part of the formula
  tittle = "Fence-line vs Truck weaning each collection day: Genus-level pairwise comparison"
  da = "All"
  vari = "weaning"
  relevel = c("Pre_weaning", "At_weaning", "Post_weaning")
# PLEASE NOTE that for geom_text_repel function in the code, you can chose which data points can have a name in the graph. You can either do: 
  # data=subset(tmp, Expr5 == 1) # will label only taxa with average exp > 5
  # data=subset(tmp, adj.P.Val < 0.05) # will label only taxa with significant adj. p-value
  
```

GENUS-LEVEL Plots
```{r}
  # Subsetting by categorical predictor (variable) and getting rid of NAs and not_match
  tmp <- subset(Log_Gen, Log_Gen$variable == vari)
  tmp <- subset(tmp, tmp$day != "All")
  tmp <- filter(tmp, taxa != "NA", taxa !="no_match")
#median(log2.phy$AveExpr)
  # re-leveling the levels of categorical variable
  tmp$day = factor(tmp$day, levels=relevel) # activate it just in case
  # creating a new column (variable) for significant p-value
  tmp$Significant <- ifelse(tmp$adj.P.Val < 0.05, "AdjP < 0.05", "Not Significant")
  
# OTHER WAY TO FILTER DATA
#Log_Phy$EffectSize <- ifelse(Log_Phy$AveExpr < 4.2057385, "ES<1", "EF>1")
#topTable_contrasts_df['taxa'] <- row.names(topTable_contrasts_df)
#Log_Phy_fil <- subset(Log_Phy,  percentile_fifty_cut%in% c("1"))       ## taxa av.exp >50th percentile
#Log_Phy_fil
#size_AveEx <- sqrt(Log_Phy_fil$AveExpr/pi)    ## size of average expression/effectsize
  
  name <- ggplot(tmp, aes(logFC, y=-log10(adj.P.Val), fill=Significant))+
  geom_point(aes(size=AveExpr), alpha=1, shape=21)+
  scale_fill_manual(values=c( "red","grey87"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  ggtitle(tittle)+         # change the title based on contrast 1,2..
  facet_wrap(~ day)+ # you can activate this if you want graphs facetted 
  geom_text_repel(data=subset(tmp, adj.P.Val < 0.05), aes(label=taxa), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))
```

Recover the name of figure
```{r}
# change the left part of the formula
  F3C <- name 
```
