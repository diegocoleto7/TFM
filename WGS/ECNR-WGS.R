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
################## IMPORT AND CLEANUP ###############
#####################################################

# Import data in BIOM format
physeq_raw <- import_biom("Celiac_noresp.biom")

# Load metadata
metadata <- read.csv("DATA_plus.cluster.csv", header = TRUE, row.names = "ID")
metadata <- sample_data(metadata)

# Format taxonomy table
physeq_raw@tax_table@.Data <- substring(physeq_raw@tax_table@.Data, 4)
colnames(physeq_raw@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Rename sample names
sample_names(physeq_raw) <- gsub("^D", "", sample_names(physeq_raw))

# Create a phyloseq object 
physeq_clean <- merge_phyloseq(physeq_raw, metadata)

# Aggregate at different taxonomic levels
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
###################### TREE CREATION ################
#####################################################

# Create phylogenetic trees
tree_species <- rtree(ntaxa(physeq_species), rooted = TRUE, tip.label = taxa_names(physeq_species))
tree_genus <- rtree(ntaxa(physeq_genus), rooted = TRUE, tip.label = taxa_names(physeq_genus))
tree_phylum <- rtree(ntaxa(physeq_phylum), rooted = TRUE, tip.label = taxa_names(physeq_phylum))

# Add trees to phyloseq objects
physeq_species <- merge_phyloseq(physeq_species, tree_species)
physeq_genus <- merge_phyloseq(physeq_genus, tree_genus)
physeq_phylum <- merge_phyloseq(physeq_phylum, tree_phylum)

#####################################################
#################### SAVE OBJECTS ###################
#####################################################

# Save phyloseq objects
saveRDS(physeq_clean, file = "C:/Users/Diego/Desktop/TFM/WGS-Celiacos/WGS-Celiacos-Noresp/WGS_species_ECNR.rds")
saveRDS(physeq_genus, file = "C:/Users/Diego/Desktop/TFM/WGS-Celiacos/WGS-Celiacos-Noresp/WGS_genus_ECNR.rds")
saveRDS(physeq_phylum, file = "C:/Users/Diego/Desktop/TFM/WGS-Celiacos/WGS-Celiacos-Noresp/WGS_phylum_ECNR.rds")

#####################################################
###################### ABUNDANCES ###################
#####################################################

# Calculate sums and select the top 4 most abundant phyla
top_phylum_counts <- taxa_sums(physeq_phylum)
top_phylum <- sort(top_phylum_counts, decreasing = TRUE)[1:4]

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
  scale_fill_carto_d(palette = "Earth") + 
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

# Calculate sums and select the top 10 most abundant genera
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
  scale_fill_carto_d(palette = "Earth") + 
  theme_minimal() +
  labs(title = "Genus", 
       x = "Genus",
       y = "log10 Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", family = "sans", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
print(plot_genus_abundance)

# Calculate sums and select the top 10 most abundant species
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
  scale_fill_carto_d(palette = "Earth") + 
  theme_minimal() +
  labs(title = "Species", 
       x = "Species",
       y = "log10 Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", family = "sans", hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
print(plot_species_abundance)

#####################################################
################## ALPHA DIVERSITY ##################
#####################################################

# Calculate alpha diversity indices
alpha_div_species <- estimate_richness(physeq_species, measures = c("Observed", "Shannon", "Simpson"))

# Assign row names based on metadata
rownames(alpha_div_species) <- rownames(metadata)

# Add cluster information to alpha diversity table
alpha_div_species_cluster <- alpha_div_species
alpha_div_species_cluster$Cluster <- metadata$Cluster

# Convert cluster data to long format for visualization
alpha_div_long_cluster <- alpha_div_species_cluster %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(cols = -c(Sample, Cluster), names_to = "Index", values_to = "Diversity")

# Perform statistical tests to evaluate alpha diversity
alpha_stats <- alpha_div_long_cluster %>%
  group_by(Index) %>%
  do({
    data <- .
    cluster1 <- data %>% filter(Cluster == 1) %>% pull(Diversity)
    cluster2 <- data %>% filter(Cluster == 2) %>% pull(Diversity)
    shapiro1 <- shapiro.test(cluster1)
    shapiro2 <- shapiro.test(cluster2)
    if (shapiro1$p.value > 0.05 & shapiro2$p.value > 0.05) {
      test <- t.test(Diversity ~ Cluster, data = data)
    } else {
      test <- wilcox.test(Diversity ~ Cluster, data = data, exact = FALSE)
    }
    data.frame(
      Index = unique(data$Index),
      Method = ifelse(shapiro1$p.value > 0.05 & shapiro2$p.value > 0.05, "t-test", "wilcox.test"),
      p.value = test$p.value
    )
  })

# View statistical test results
View(alpha_stats)

# Create alpha diversity boxplot by Cluster
plot_alpha_species <- ggplot(alpha_div_long_cluster, aes(x = as.factor(Cluster), y = Diversity, fill = as.factor(Cluster))) +
  geom_boxplot() +
  facet_wrap(~ Index, scales = "free") +
  labs(title = "WGS",
       x = "Cluster",
       y = "Diversity",
       fill = "Cluster") +
  scale_fill_manual(values = c("1" = "#A16928", "2" = "#2887A1")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14))
print(plot_alpha_species)

# Calculate mean alpha diversity for each index by cluster
alpha_mean_values <- aggregate(Diversity ~ Cluster + Index, data = alpha_div_long_cluster, FUN = mean)

# Print mean values
print(alpha_mean_values)

#####################################################
################## BETA DIVERSITY ###################
#####################################################

# Transform data to compositional
physeq_compositional <- microbiome::transform(physeq_species, "compositional")
physeq_compositional@sam_data$Cluster <- as.factor(physeq_compositional@sam_data$Cluster)

# Calculate distances
dist_bray <- phyloseq::distance(physeq_compositional, method = "bray")
dist_unifrac <- phyloseq::distance(physeq_compositional, method = "unifrac", weighted = TRUE)

# Perform PERMANOVA
permanova_bray <- adonis2(dist_bray ~ Cluster, data = as(sample_data(physeq_compositional), "data.frame"))
permanova_unifrac <- adonis2(dist_unifrac ~ Cluster, data = as(sample_data(physeq_compositional), "data.frame"))

# Print results
print("PERMANOVA Bray-Curtis results")
print(permanova_bray)
print("PERMANOVA UniFrac results")
print(permanova_unifrac)

# UniFrac PCoA plot
ordinated_taxa_unifrac <- ordinate(physeq = physeq_compositional, method = "PCoA", distance = "unifrac", weighted = TRUE)
p_unifrac_pcoa <- plot_ordination(physeq_compositional, ordinated_taxa_unifrac, color = "Cluster", title = "WGS ") +
  geom_point(size = 4) +
  theme_bw() +
  scale_color_manual(values = c("1" = "#A16928", "2" = "#2887A1")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
print(p_unifrac_pcoa)

#####################################################
############### ABUNDANCE WITH DESEQ2 ###############
#####################################################

# Correct OTU names
tax_df <- as.data.frame(tax_table(physeq_species))
new_taxa_names <- paste(tax_df$Genus, tax_df$Species, sep = "_")
taxa_names(physeq_species) <- new_taxa_names

# Filter taxa
flist <- filterfun(kOverA(k = 3, A = 1))
data.fill <- filter_taxa(physeq_species, flist)
data_f <- prune_taxa(data.fill, physeq_species)

# Prepare data for DESeq2
sample_data(data_f)$Cluster <- as.factor(sample_data(data_f)$Cluster)
dds <- phyloseq_to_deseq2(data_f, ~ Cluster)

# Run differential abundance analysis
dds <- DESeq(dds)
res <- results(dds)
res_ordered <- res[order(res$padj), ]

# Show main results
head(res_ordered)
