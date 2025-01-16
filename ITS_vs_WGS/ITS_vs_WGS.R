library(phyloseq)
library(ggplot2)
library(microbiome)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggrepel)
library(VennDiagram)
library(grid)

# Load the WGS phyloseq object at the Genus level
ps_shotgun <- readRDS("path/to/genus_physeq_WGS.rds")

# Load the ITS phyloseq object at the Genus level
ps_its <- readRDS("path/to/genus_physeq_ITS.rds")

# Extract column names
otu_cols_ps_its <- colnames(otu_table(ps_its))
otu_cols_ps_shotgun <- colnames(otu_table(ps_shotgun))

# Identify sample names present in both and prune
in_both <- intersect(otu_cols_ps_its, otu_cols_ps_shotgun)
ps_shotgun <- prune_samples(sample_names(ps_shotgun) %in% in_both, ps_shotgun)
ps_its <- prune_samples(sample_names(ps_its) %in% in_both, ps_its)

# Remove  zeros
ps_shotgun <- prune_taxa(taxa_sums(ps_shotgun) > 0, ps_shotgun)
ps_its <- prune_taxa(taxa_sums(ps_its) > 0, ps_its)

# Load the shotgun phyloseq object at the Species level
especies_shotgun <- readRDS("path/to/species_physeq_WGS.rds")

# Load the ITS phyloseq object at the Species level
especies_its <- readRDS("path/to/species_physeq_ITS.rds")

# Prune to match the samples
especies_shotgun <- prune_samples(sample_names(especies_shotgun) %in% sample_names(ps_shotgun), especies_shotgun)
especies_its <- prune_samples(sample_names(especies_its) %in% sample_names(ps_its), especies_its)


######################################################################

# Identify basic statistics for ITS: max, min, mean, total reads
max_its <- max(sample_sums(ps_its))
min_its <- min(sample_sums(ps_its))
mean_its <- mean(sample_sums(ps_its))
sum_its <- sum(sample_sums(ps_its))

cat("ITS - Max reads:\n", max_its, "\n")
cat("ITS - Min reads:\n", min_its, "\n")
cat("ITS - Mean reads:\n", mean_its, "\n")
cat("ITS - Total reads:\n", sum_its, "\n")

# Identify basic statistics for Shotgun: max, min, mean, total reads
max_shotgun <- max(sample_sums(ps_shotgun))
min_shotgun <- min(sample_sums(ps_shotgun))
mean_shotgun <- mean(sample_sums(ps_shotgun))
sum_shotgun_val <- sum(sample_sums(ps_shotgun))

cat("Shotgun - Max reads:\n", max_shotgun, "\n")
cat("Shotgun - Min reads:\n", min_shotgun, "\n")
cat("Shotgun - Mean reads:\n", mean_shotgun, "\n")
cat("Shotgun - Total reads:\n", sum_shotgun_val, "\n")

# Show total reads in the species-level Shotgun object
cat("Shotgun (Species object) - Total reads:\n", sum(sample_sums(especies_shotgun)), "\n")

# Calculate the percentage of ITS reads assigned to species level
per_its <- (100 * sum(sample_sums(especies_its))) / sum_its
cat("Percentage of ITS reads assigned to species:\n", per_its, "\n")

# Calculate the percentage of Shotgun reads assigned to species level
per_shotgun <- (100 * sum(sample_sums(especies_shotgun))) / sum_shotgun_val
cat("Percentage of Shotgun reads assigned to species:\n", per_shotgun, "\n")

############################################################################
# Taxonomic agglomeration  and trasform
its_phyl_AI <- tax_glom(ps_its, taxrank = "Phylum")
ps_its_compositional <- microbiome::transform(its_phyl_AI, transform = "compositional")

# Calculate average abundance at the Phylum level for ITS
its_data <- psmelt(ps_its_compositional) %>%
  group_by(Phylum) %>%
  summarise(Abundance = mean(Abundance)) %>%
  mutate(Abundance = Abundance * 100)

# Taxonomic agglomeration and transform
shotgun_phyl_AI <- tax_glom(ps_shotgun, taxrank = "Phylum")
ps_shotgun_compositional <- microbiome::transform(shotgun_phyl_AI, transform = "compositional")

# Calculate average abundance at the Phylum level for Shotgun
shotgun_data <- psmelt(ps_shotgun_compositional) %>%
  group_by(Phylum) %>%
  summarise(Abundance = mean(Abundance)) %>%
  mutate(Abundance = Abundance * 100)

# Arrange phyla in a desired order for ITS
its_data$Phylum <- factor(
  its_data$Phylum, 
  levels = c("Ascomycota", "Chytridiomycota", "Mortierellomycota", "Mucoromycota", 
             "Sanchytriomycota", "Basidiomycota")
)

# Create pie chart for ITS
plot_its_pie <- ggplot(its_data, aes(x = "", y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  coord_polar("y") +
  ggtitle("Phylum - ITS") +
  scale_fill_carto_d(palette = "Earth") + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  ) +
  geom_text(
    aes(label = ifelse(Abundance > 2, paste0(round(Abundance, 1), "%"), "")), 
    position = position_stack(vjust = 0.5), 
    size = 5
  )

plot_its_pie

# Arrange phyla in a desired order for Shotgun
shotgun_data$Phylum <- factor(
  shotgun_data$Phylum, 
  levels = c("Ascomycota", "Microsporidia", "Mucoromycota", "Basidiomycota")
)

# Create pie chart for Shotgun
plot_shotgun_pie <- ggplot(shotgun_data, aes(x = "", y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.2) +
  coord_polar("y") +
  ggtitle("Phylum - WGS") +
  scale_fill_carto_d(palette = "Earth") + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  ) +
  geom_text(
    aes(label = ifelse(Abundance > 2, paste0(round(Abundance, 1), "%"), "")), 
    position = position_stack(vjust = 0.5), 
    size = 5
  )

plot_shotgun_pie

#########################################################################
# Identify genera in ITS
genus_its <- unique(tax_table(ps_its)[, "Genus"])

# Identify genera in Shotgun
genus_shotgun <- unique(tax_table(ps_shotgun)[, "Genus"])

# Create list for Venn diagram
venn_list <- list(WGS = genus_shotgun, ITS = genus_its)

# Generate the Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("WGS", "ITS"),
  filename = NULL,
  fill = c("#A16928", "#2887A1"),
  alpha = 0.8,
  lty = "solid",
  lwd = 0.5,
  cex = 1.6,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 1.2,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-30, 30),
  cat.dist = c(0.05, 0.05),
  main = "Diagrama de Venn - Géneros",
  main.cex = 2,
  main.fontface = "bold",
  main.fontfamily = "sans",
  sub = "Celiacos",
  sub.cex = 1.5,
  sub.fontface = "bold",
  sub.fontfamily = "sans"
)

grid.draw(venn.plot)

##########################################################

# Identify common genera between ITS and Shotgun
common_genera <- intersect(genus_shotgun, genus_its)

# Convert both objects to compositional data
percentaje_its <- transform(ps_its, "compositional")
percentaje_shotgun <- transform(ps_shotgun, "compositional")

# Extract  abundances for ITS
otu_its_df <- as.data.frame(taxa_sums(percentaje_its))
tax_its_df <- as.data.frame(tax_table(percentaje_its))
otu_its_df$OTU_ID <- rownames(otu_its_df)
tax_its_df$OTU_ID <- rownames(tax_its_df)
merged_its_df <- merge(otu_its_df, tax_its_df, by = "OTU_ID")
genus_abundance_its <- merged_its_df %>%
  group_by(Genus) %>%
  summarise(ITS_Abundance = `taxa_sums(percentaje_its)`)

# Extract  abundances for Shotgun
otu_shotgun_df <- as.data.frame(taxa_sums(percentaje_shotgun))
tax_shotgun_df <- as.data.frame(tax_table(percentaje_shotgun))
otu_shotgun_df$OTU_ID <- rownames(otu_shotgun_df)
tax_shotgun_df$OTU_ID <- rownames(tax_shotgun_df)
merged_shotgun_df <- merge(otu_shotgun_df, tax_shotgun_df, by = "OTU_ID")
genus_abundance_shotgun <- merged_shotgun_df %>%
  group_by(Genus) %>%
  summarise(WGS_Abundance = `taxa_sums(percentaje_shotgun)`)

# Filter to keep only common genera in both datasets
genus_abundance_its_filtered <- genus_abundance_its %>%
  filter(Genus %in% common_genera)
genus_abundance_shotgun_filtered <- genus_abundance_shotgun %>%
  filter(Genus %in% common_genera)

# Merge both dataframes by Genus
combined_data <- merge(genus_abundance_its_filtered, genus_abundance_shotgun_filtered, by = "Genus")

# Perform Pearson correlation test
cor_result <- cor.test(combined_data$ITS_Abundance, combined_data$WGS_Abundance, method = "pearson")
print(cor_result)

# Prepare correlation label
corr_label <- paste0(
  "r = ", round(cor_result$estimate, 2), 
  "\np-value = ", signif(cor_result$p.value, 2)
)

# Plot scatter diagram with correlation
ggplot(combined_data, aes(x = ITS_Abundance, y = WGS_Abundance, label = Genus)) +
  geom_point(size = 3, alpha = 0.8, aes(color = Genus)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_text_repel(size = 3) +
  scale_x_sqrt() +
  scale_y_sqrt() +
  labs(
    title = "Correlación de Abundancia entre Métodos",
    x = "Abundancia Relativa (ITS)",
    y = "Abundancia Relativa (Shotgun)",
    color = "Género"
  ) +
  annotate("text", x = Inf, y = Inf, label = corr_label, hjust = 2.5, vjust = 2.5, size = 3) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
