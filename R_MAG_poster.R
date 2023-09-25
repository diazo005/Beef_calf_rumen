## MAGs poster
## It will present the prelimnary results of 697 No-DeReplicated MAGs obtained with 50% completeness and 10% contamination

## load packages
source("Packages.R")

## Import data
reassembly <- read_excel("MAG_Master.xlsx", sheet = "reassembly_log")
reassembly <- as.data.frame(reassembly)
reassembly <- filter(reassembly, reassembly$sample_type == "Sample")
reassembly$collection_day = factor(reassembly$collection_day, levels=c("Pre_weaning","At_weaning","Post_weaning"))
levels(reassembly$collection_day) <- c('Pre weaning', 'Weaning day', 'Post weaning')

## Create scatter plot with marginal histograms/desity plots
p <- ggscatterhist(
  reassembly, x = "contamination", y = "completeness",
  color = "collection_day", size = 3, alpha = 0.6,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.params = list(fill = "collection_day", color = "black", size = 0.2)
)     

## statistics of completeness and contamination
describeBy(reassembly$completeness)
describeBy(reassembly$contamination)


## Import data
taxa <- read_excel("MAG_Master.xlsx", sheet = "gtdbtk_log")
taxa <- as.data.frame(taxa)
taxa <- filter(taxa, taxa$sample_type == "Sample")
taxa$collection_day = factor(taxa$collection_day, levels=c("Pre_weaning","At_weaning","Post_weaning"))
levels(taxa$collection_day) <- c('Pre weaning', 'Weaning day', 'Post weaning')

## Create abundance barplot with counts (not percentage)
p2 <- ggplot(taxa %>% count(collection_day, Phylum) %>%    # Group by region and species, then count number in each group
         mutate(pct=n/sum(n),               # Calculate percent within each region
                ypos = cumsum(n) - 0.5*n),  # Calculate label positions
       aes(collection_day, n, fill=Phylum)) +
  geom_bar(stat="identity") +
  ylab("MAGs (counts)")+theme_bw()+
  ggthemes::scale_fill_tableau("Tableau 20")+
  ggtitle("Number of MAGs per collection day classified at phylum level")         # change the title based on contrast 1,2..
  
## Import data
samples <- read_excel("MAG_Master.xlsx", sheet = "Samples")
samples <- as.data.frame(samples)

## Create table with summary statistics of MAG building process (Meas+SD)
describeBy(samples$reads_raw,samples$collection_day)
describeBy(samples$reads_nonhost,samples$collection_day)
describeBy(samples$contigs_n,samples$collection_day)
describeBy(samples$contigs_avg_bp,samples$collection_day)
describeBy(samples$`bins_pre-refinement_avg`,samples$collection_day)

