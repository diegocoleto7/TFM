library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ape)
library(microbiome)
library(rcartocolor)
library(DESeq2)
library(genefilter)

#####################################################
############## IMPORT AND CLEANUP ###################
#####################################################

# Define paths
feature_table_path <- "C:/Users/Diego/Desktop/TFM/ITS-AI/nochimera-AI/exported-feature-table/feature-table.tsv"
taxonomy_path <- "C:/Users/Diego/Desktop/TFM/ITS-AI/nochimera-AI/exported-taxonomy/taxonomy.tsv"

# Import feature table
otu_table <- read.table(feature_table_path, header = TRUE, row.names = 1, sep = "\t", skip = 1, check.names = FALSE, comment.char = "")
otu_table <- as.matrix(otu_table)
OTU <- otu_table(otu_table, taxa_are_rows = TRUE)

# Import taxonomy
taxonomy <- read.table(taxonomy_path, header = TRUE, row.names = 1, sep = "\t")
taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";", fill = "right", extra = "drop")
tax_matrix <- as.matrix(taxonomy)
TAX <- tax_table(tax_matrix)

# Create the phyloseq object
physeq <- phyloseq(OTU, TAX)

# Filter by confidence
physeq_clean <- prune_taxa(rownames(taxonomy[taxonomy$Confidence >= 0.95, ]), physeq)

# Clean taxonomy
physeq_clean@tax_table@.Data <- substring(physeq_clean@tax_table@.Data, 4)

# Aggregate by different taxonomic levels
physeq_phylum <- tax_glom(physeq_clean, taxrank = "Phylum")
physeq_genus <- tax_glom(physeq_clean, taxrank = "Genus")
physeq_species <- tax_glom(physeq_clean, taxrank = "Species")

# Remove OTUs and samples with zero total counts
physeq_phylum <- prune_taxa(taxa_sums(physeq_phylum) > 0, physeq_phylum)
physeq_phylum <- prune_samples(sample_sums(physeq_phylum) > 0, physeq_phylum)

physeq_genus <- prune_taxa(taxa_sums(physeq_genus) > 0, physeq_genus)
physeq_genus <- prune_samples(sample_sums(physeq_genus) > 0, physeq_genus)

physeq_species <- prune_taxa(taxa_sums(physeq_species) > 0, physeq_species)
physeq_species <- prune_samples(sample_sums(physeq_species) > 0, physeq_species)

#####################################################
################### METADATA ########################
#####################################################

# Create metadata based on sample names for species
sample_names_data <- sample_names(physeq_species)
Visit <- ifelse(grepl("^V1_", sample_names_data), "Visit1",
                ifelse(grepl("^V3_", sample_names_data), "Visit2", NA))
ID <- gsub("^V[13]_", "", sample_names_data)
metadata <- data.frame(row.names = sample_names_data, Visit = Visit, ID = ID)

metadata <- sample_data(metadata)
physeq_species <- merge_phyloseq(physeq_species, metadata)

# Create metadata based on sample names for genus
sample_names_genus <- sample_names(physeq_genus)
Visit_genus <- ifelse(grepl("^V1_", sample_names_genus), "Visit1",
                      ifelse(grepl("^V3_", sample_names_genus), "Visit2", NA))
ID_genus <- gsub("^V[13]_", "", sample_names_genus)
metadata_genus <- data.frame(row.names = sample_names_genus, Visit = Visit_genus, ID = ID_genus)

metadata_genus <- sample_data(metadata_genus)
physeq_genus <- merge_phyloseq(physeq_genus, metadata_genus)

# Create metadata based on sample names for phylum
sample_names_phylum <- sample_names(physeq_phylum)
Visit_phylum <- ifelse(grepl("^V1_", sample_names_phylum), "Visit1",
                       ifelse(grepl("^V3_", sample_names_phylum), "Visit2", NA))
ID_phylum <- gsub("^V[13]_", "", sample_names_phylum)
metadata_phylum <- data.frame(row.names = sample_names_phylum, Visit = Visit_phylum, ID = ID_phylum)

metadata_phylum <- sample_data(metadata_phylum)
physeq_phylum <- merge_phyloseq(physeq_phylum, metadata_phylum)

#####################################################
################### TREE CREATION ###################
#####################################################

# Create and add phylogenetic trees
tree_species <- rtree(ntaxa(physeq_species), rooted = TRUE, tip.label = taxa_names(physeq_species))
tree_genus <- rtree(ntaxa(physeq_genus), rooted = TRUE, tip.label = taxa_names(physeq_genus))
tree_phylum <- rtree(ntaxa(physeq_phylum), rooted = TRUE, tip.label = taxa_names(physeq_phylum))

physeq_species <- merge_phyloseq(physeq_species, tree_species)
physeq_genus <- merge_phyloseq(physeq_genus, tree_genus)
physeq_phylum <- merge_phyloseq(physeq_phylum, tree_phylum)

# Save phyloseq objects
saveRDS(physeq_species, file = "C:/Users/Diego/Desktop/TFM/ITS-AI/ITS_species_AI.rds")
saveRDS(physeq_genus, file = "C:/Users/Diego/Desktop/TFM/ITS-AI/ITS_genus_AI.rds")
saveRDS(physeq_phylum, file = "C:/Users/Diego/Desktop/TFM/ITS-AI/ITS_phylum_AI.rds")

#####################################################
################### ABUNDANCES ######################
#####################################################

# Calculate sums and select the 10 most abundant phyla
top_phylum_counts <- taxa_sums(physeq_phylum)
top_phylum <- sort(top_phylum_counts, decreasing = TRUE)[1:10]

# Get the taxonomy associated with the top abundant phyla
taxonomy_phylum_top <- tax_table(physeq_phylum)[names(top_phylum), ]

# Create a dataframe for visualization
df_phylum <- as.data.frame(taxonomy_phylum_top)
df_phylum$Abundance <- top_phylum
df_phylum$Phylum <- factor(df_phylum$Phylum, levels = df_phylum %>% 
                             group_by(Phylum) %>% 
                             summarize(TotalAbundance = sum(Abundance)) %>% 
                             arrange(desc(TotalAbundance)) %>% 
                             pull(Phylum))

# Plot phylum abundance
plot_phylum_abundance <- ggplot(df_phylum, aes(x = Phylum, y = log10(Abundance), fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_carto_d(palette = "Geyser") + 
  theme_minimal() +
  labs(title = "Phylum", 
       x = "Phylum",
       y = "log10 Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", family = "sans", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
print(plot_phylum_abundance)

# Calculate sums and select the 10 most abundant genera
top_genus_counts <- taxa_sums(physeq_genus)
top_genus <- sort(top_genus_counts, decreasing = TRUE)[1:10]

# Get the taxonomy associated with the top abundant genera
taxonomy_genus_top <- tax_table(physeq_genus)[names(top_genus), ]

# Create a dataframe for visualization
df_genus <- as.data.frame(taxonomy_genus_top)
df_genus$Abundance <- top_genus
df_genus$Genus <- factor(df_genus$Genus, levels = df_genus %>% 
                           group_by(Genus) %>% 
                           summarize(TotalAbundance = sum(Abundance)) %>% 
                           arrange(desc(TotalAbundance)) %>% 
                           pull(Genus))

# Plot genus abundance
plot_genus_abundance <- ggplot(df_genus, aes(x = Genus, y = log10(Abundance), fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_carto_d(palette = "Geyser") + 
  theme_minimal() +
  labs(title = "Género", 
       x = "Género",
       y = "log10 Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", family = "sans", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
print(plot_genus_abundance)

# Calculate sums and select the 10 most abundant species
top_species_counts <- taxa_sums(physeq_species)
top_species <- sort(top_species_counts, decreasing = TRUE)[1:10]

# Get the taxonomy associated with the top abundant species
taxonomy_species_top <- tax_table(physeq_species)[names(top_species), ]

# Create a dataframe for visualization
df_species <- as.data.frame(taxonomy_species_top)
df_species$Abundance <- top_species
df_species$Species <- factor(df_species$Species, levels = df_species %>% 
                               group_by(Species) %>% 
                               summarize(TotalAbundance = sum(Abundance)) %>% 
                               arrange(desc(TotalAbundance)) %>% 
                               pull(Species))

# Plot species abundance
plot_species_abundance <- ggplot(df_species, aes(x = Species, y = log10(Abundance), fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_carto_d(palette = "Geyser") + 
  theme_minimal() +
  labs(title = "Especie", 
       x = "Especie",
       y = "log10 Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", family = "sans", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
print(plot_species_abundance)

#####################################################
################### ALPHA DIVERSITY #################
#####################################################

# Calculate alpha diversity indices
alpha_div_species <- estimate_richness(physeq_species, measures = c("Observed", "Shannon", "Simpson"))
rownames(alpha_div_species) <- rownames(metadata)

# Add visit information to the alpha diversity table
alpha_div_species_visit <- alpha_div_species
alpha_div_species_visit$Visit <- metadata$Visit

# Table for plotting
alpha_div_long_visit <- alpha_div_species_visit %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(cols = -c(Sample, Visit), names_to = "Index", values_to = "Diversity")

# Perform statistical tests to evaluate alpha diversity
alpha_stats <- alpha_div_long_visit %>%
  group_by(Index) %>%
  do({
    data <- .
    visit1 = data %>% filter(Visit == "Visit1") %>% pull(Diversity)
    visit2 = data %>% filter(Visit == "Visit2") %>% pull(Diversity)
    shapiro1 <- shapiro.test(visit1)
    shapiro2 <- shapiro.test(visit2)
    
    if (shapiro1$p.value > 0.05 & shapiro2$p.value > 0.05) {
      test <- t.test(Diversity ~ Visit, data = data) 
    } else {
      test <- wilcox.test(Diversity ~ Visit, data = data, exact = FALSE) 
    }
    data.frame(
      Index = unique(data$Index),
      Method = ifelse(shapiro1$p.value > 0.05 & shapiro2$p.value > 0.05, "t-test", "wilcox.test"),
      p.value = test$p.value
    )
  })

# Create alpha diversity plot by Visit
p_alpha_cluster <- ggplot(alpha_div_long_visit, aes(x = as.factor(Visit), y = Diversity, fill = as.factor(Visit))) +
  geom_boxplot() +
  facet_wrap(~ Index, scales = "free") +
  labs(title = "ITS",
       x = "Visit",
       y = "Diversity" ,
       fill = "Visit") +
  scale_fill_manual(values = c("Visit1" = "#008080", "Visit2" = "#CA562C")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14)) 
p_alpha_cluster

# Calculate the mean of each index by cluster
mean_values <- aggregate(Diversity ~ Visit + Index, data = alpha_div_long_visit, FUN = mean)

# Print results
print(mean_values)

#####################################################
################### BETA DIVERSITY ##################
#####################################################

# Transform to compositional data and calculate distances
physeq_compositional <- microbiome::transform(physeq_species, "compositional")
physeq_compositional@sam_data$Visit <- as.factor(physeq_compositional@sam_data$Visit)

dist_bray <- phyloseq::distance(physeq_compositional, method = "bray")
dist_unifrac <- phyloseq::distance(physeq_compositional, method = "unifrac", weighted = TRUE)

# PERMANOVA
permanova_bray <- adonis2(dist_bray ~ Visit, data = as(sample_data(physeq_compositional), "data.frame"))
permanova_unifrac <- adonis2(dist_unifrac ~ Visit, data = as(sample_data(physeq_compositional), "data.frame"))

print("PERMANOVA Bray-Curtis:")
print(permanova_bray)
print("PERMANOVA UniFrac:")
print(permanova_unifrac)

# UniFrac PCoA plot
ordinated_taxa_unifrac <- ordinate(physeq = physeq_compositional, method = "PCoA", distance = "unifrac", weighted = TRUE)
p_unifrac_pcoa <- plot_ordination(physeq_compositional, ordinated_taxa_unifrac, color = "Visit", title = "ITS") +
  geom_point(size = 2) +
  theme_bw() +
  scale_color_manual(values = c("Visit1" = "#008080", "Visit2" = "#D7784A")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
print(p_unifrac_pcoa)

#####################################################
################## DESEQ2 ABUNDANCES ################
#####################################################

taxa_names(physeq_species) <- tax_table(physeq_species)[,"Species"]

flist <- filterfun(kOverA(k = 3, A = 10))
data.fill <- filter_taxa(physeq_species, flist)
data_f <- prune_taxa(data.fill, physeq_species)

# Convert the phyloseq object to a DESeq2 object
sample_data(data_f)$Visit <- as.factor(sample_data(data_f)$Visit)

cont <- otu_table(data_f)
cont <- cont + 1
otu_table(data_f) <- cont
dds <- phyloseq_to_deseq2(data_f, ~ Visit)

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# See the top results
head(res_ordered)
