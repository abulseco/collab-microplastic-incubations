# MICROPLASTICS PHYLOSEQ ANALYSIS
# SEQUENCES ANALYZED IN DADA2 (microplastics-DADA2.R)
# LAST UPDATED: 12/30/25 BY ANB
# Code revised from "microplastics-PHYLOSEQ-UPDATE.R" 

## LIBRARIES-------------------------------------------
library(phyloseq); library(data.table); library(ggplot2); library(dplyr)
library(viridis); library(RColorBrewer); library(indicspecies)

## FORMATTING FUNCTIONS--------------------------------
pretty.theme <- function(){
  theme_bw() +
    theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black"),
          axis.text.y=element_text(size = 18, color="black"),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          text=element_text(size=18),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA),
          axis.line = element_line(colour = "black"))
}

# Loading into Phyloseq----
mat = read.table("INPUT-FILES/OTU-sample-IDS.txt", header = TRUE, sep = "\t", row.names = 1)
tax = read.table("INPUT-FILES/ASVs_taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1)
meta = read.csv("INPUT-FILES/sample-metadata.csv", header = TRUE)
# meta = read.table("sample-metadata.txt", header = TRUE, sep = "\t", row.names = 1)
# tree = read_tree("asv-famsa-trimal.tre") # Still need to make the tree

mat = as.matrix(mat)
tax = as.matrix(tax)

OTU = otu_table(mat, taxa_are_rows = TRUE)
TAX = tax_table(tax)
META = sample_data(meta)
# sample_names(OTU)
# sample_names(TAX)
# tail(META)
# head(TAX)
# head(OTU)

phy = phyloseq(OTU,TAX,META)
phy
sample_names(phy)
# NOTE: Eventually need to fix this as the sample IDs are still "sa1" 

# Learn about phyloseq object
ntaxa(phy) # 5107
nsamples(phy) # Check that it is what you expect - 51
sample_variables(phy) # Should match your metadata file
taxa_names(phy)
sample_names(phy)
sample_sums(phy)

# Look at statistics related to your data
num.reads = sum(sample_sums(phy))
lowest.sam = sort(sample_sums(phy)) 
mean.seq = mean(sample_sums(phy))  
std.seq = sd(sample_sums(phy))/sqrt(51) 
med.seq = median(sample_sums(phy))
phy.summary <- data.frame(num.reads, mean.seq, std.seq, med.seq)
phy.summary     

seq.dt = data.table(as(sample_data(phy), "data.frame"),
                    TotalReads = sample_sums(phy), keep.rownames = TRUE)
seq.dt # information per sample after filtering out taxa 
View(seq.dt)

# FILTERING----
# Get rid of Mitochondria  
phy %>%
  subset_taxa(Family != c("Mitochondria", "Chloroplast")) -> phy.f
phy # BEFORE filtering: 5107 taxa, 51 samples
phy.f # AFTER filtering: 2894 taxa, 51 samples

# If you want a version that keeps chloroplasts
phy %>%
  subset_taxa(Family != "Mitochondria") -> phy.f.wChloro
phy # BEFORE filtering: 5107 taxa, 51 samples
phy.f.wChloro # AFTER filtering: 2894 taxa, 51 samples

# remove blank sample
phy.f <- subset_samples(phy.f, Sample_Type != "Blank")
phy.f # Should be 50 samples if the blank was removed

# For microplastics only
# Need to do this in two steps, ignore the "1" moving forward
phyf.micro1 <- phy.f %>% 
  subset_samples(Sample_Type != "Aggregate")

phyf.micro <- phyf.micro1 %>% 
  subset_samples(Sample_Type != "Natural Particles")

# Subset the data to get rid of aggregates
phyf.noagg <- phy.f %>%
  subset_samples(Sample_Type != "Aggregate")
sample_data(phyf.noagg)

# Subset the data to look at natural particles only
phyf.noagg.nat <- phyf.micro1 %>% 
  subset_samples(Sample_Type == "Natural Particles")
sample_data(phyf.noagg.nat)

# Filter out low abundances
lowabundnames = filter_taxa(phy.f, function(x) sum(x) > 1)
phy.f.nolow = prune_taxa(lowabundnames, phy.f) 
ntaxa(phy.f.nolow) # 2847 remain

# Calculate relative abundance
per.f = transform_sample_counts(phy.f, function (x) x/sum(x)*100)
per.f.nolow = transform_sample_counts(phy.f.nolow, function (x) x/sum(x)*100)

# Removing potentially contaminated sample from the dataset
# Define the sample(s) you want to remove
# samples_to_remove <- "sa9"  # Replace with your specific sample name(s)
# 
# # Filter out the sample(s)
# phyf.noagg <- prune_samples(!(sample_names(phyf.noagg) %in% samples_to_remove), phyf.noagg)

# Assign order
# # Assign order
desired_order_days <- c("2", "5", "9", "12", "16", "23")  # Replace with your actual day names
desired_order_sample_types <- c("Aggregate", "Natural Particles", "Microplastic ")

# # Reorder Day
sample_data(phy.f)$Day <- factor(sample_data(phy.f)$Day, levels = desired_order_days)
sample_data(phy.f)$Sample_Type <- factor(sample_data(phy.f)$Sample_Type, levels = desired_order_sample_types)

# Plots----
## Ordination no aggregates----
BC_distance <- phyloseq::distance(phy.f, "bray")
bcOrd <- ordinate(phy.f, "PCoA", BC_distance)
# bcOrd <- ordinate(phy.f, "NMDS", BC_distance)
# plot_scree(bcOrd)

p1 <- plot_ordination(phy.f, bcOrd, color = "Day", shape = "Sample_Type") +
  geom_point(size = 8, alpha = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  # geom_text(aes(label = Sample_ID), color = "black") + # Identify outlier sample to reseq 
  pretty.theme() +
  labs(shape = "Sample Type")
p1

## Barplot----
## Barplot averaged by sample (to fix the barplot)---------
# phyf.noagg.rel <- transform_sample_counts(phyf.noagg, function(x) (x / sum(x))*100)

# Sum the abundances for each ASV across all samples
asv_totals <- taxa_sums(phy.f)

# Get the names of the top 20 most abundant ASVs
top_20_asvs <- names(sort(asv_totals, decreasing = TRUE)[1:20])

# Subset the phyloseq object to include only the top 20 ASVs
physeq_top20 <- prune_taxa(top_20_asvs, phy.f)

# Then do we need to calculate relative abundance? 
phy.f.top20 <- transform_sample_counts(physeq_top20, function(x) (x / sum(x))*100)

# Melt the phyloseq object to long format
physeq_melt <- psmelt(phy.f.top20)

# Calculate mean abundance for each ASV across replicates
physeq_mean <- physeq_melt %>%
  group_by(Day, Family, Sample_Type) %>%
  summarize(Mean_Abundance = mean(Abundance))

# Check the summarized data
print(physeq_mean)

# Change the order
desired_order_days <- c("2", "5", "9", "12", "16", "23")  # Replace with your actual day names
desired_order_sample_types <- c("Aggregate", "Natural Particles", "Microplastic ")

physeq_mean$Day <- factor(physeq_mean$Day, levels = desired_order_days)
physeq_mean$Sample_Type <- factor(physeq_mean$Sample_Type, levels = desired_order_sample_types)

colors_20 <- c("#1B9E77", "#d98302", "#7570B3", "#E7298A", "#15425e", "#E6AB02", "#8e998d",
               "#c22131", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
               "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

p <- ggplot(physeq_mean, aes(x = Day, y = Mean_Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", color = "NA") +
  labs(x = "Day of Sampling", y = "Mean Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~Sample_Type) +
  scale_fill_manual(values = colors_20) +
  pretty.theme()
p

p2 <- ggplot(physeq_mean, aes(x = Day, y = Mean_Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(x = "Day of Sampling", y = "Mean Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~Sample_Type) +
  scale_fill_manual(values = colors_20) +
  pretty.theme()
p2

## Indicator analysis----
# Calculate relative abundances
phy.f.rel <- transform_sample_counts(phy.f, function(x) x / sum(x))

# Filter ASVs to include only those with at least 1% relative abundance in any sample
phy.f.filtered <- prune_taxa(taxa_sums(phy.f.rel) > 0.01, phy.f.rel)

# Extract the OTU table and sample data from the filtered phyloseq object
otu_table <- otu_table(phy.f.filtered)
sample_data <- sample_data(phy.f.filtered)

# Ensure OTU table is a matrix
otu_matrix <- as(otu_table(phy.f.filtered), "matrix")

# Ensure the sample data has the grouping variable (e.g., Sample_Type)
grouping_variable <- sample_data(phy.f.filtered)$Sample_Type

# Ensure that the OTU table has taxa as rows and samples as columns
if (taxa_are_rows(phy.f.filtered)) {
  otu_matrix <- t(otu_matrix)
}

# Check dimensions to ensure they match
print(dim(otu_matrix))
print(length(grouping_variable))

# Perform Indicator Species Analysis
isa_result <- multipatt(otu_matrix, grouping_variable, func = "IndVal.g", control = how(nperm=999))

# View the results
summary(isa_result)

# Extract significant indicator species (e.g., p < 0.05)
significant_species <- isa_result$sign[isa_result$sign$p.value < 0.05,]

# View significant species
print(significant_species)

# Extract significant OTUs
significant_otus <- rownames(significant_species)

# Rank significant species by p-value and select the top 20
significant_species <- significant_species[order(significant_species$p.value), ]
top_20_species <- significant_species[1:20, ]

# Extract significant OTUs
top_20_otus <- rownames(top_20_species)
print(top_20_otus)

# Subset the phyloseq object to include only the top 20 significant OTUs
phy_sig <- prune_taxa(top_20_otus, phy.f.filtered)
phy_sig_glom <- tax_glom(phy_sig, taxrank = "Family")

# Melt it!
df_phy_sig <- psmelt(phy_sig)

# Preparing for the plot (order and colors)
df_phy_sig$Day <- factor(df_phy_sig$Day, levels = desired_order_days)
df_phy_sig$Sample_Type <- factor(df_phy_sig$Sample_Type, levels = desired_order_sample_types)
colors_20 <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))

# Create the stacked bar plot
ggplot(df_phy_sig, aes(x = Day, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Sample", y = "Relative Abundance", fill = "Family") +
  facet_grid(~Sample_Type) +  # Facet by Sample_Type and Day
  scale_fill_manual(values = colors_20) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))




