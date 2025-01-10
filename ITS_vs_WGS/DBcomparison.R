library(phyloseq)
library(ggplot2)

# Read RDS files
genus_old <- readRDS("/path/to/genus_WGS_celiacos.rds")
genus_new <- readRDS("/path/to/Comparacion_WGS_celiacos_new.rds")

species_old <- readRDS("/path/to/phyloseq_WGS_celiacos.rds")
species_new <- readRDS("/path/to/phyloseq_WGS_celiacos_newdb.rds")

phylum_old <- readRDS("/path/to/phylum_WGS_celiacos.rds")
phylum_new <- readRDS("/path/to/phylum_WGS_celiacos_new.rds")

# Function to count the number of unique taxa at each taxonomic level
count_taxa <- function(physeq_obj) {
  taxa_table <- tax_table(physeq_obj)
  apply(taxa_table, 2, function(x) length(unique(na.omit(x))))
}

# Calculate the number of taxa for each database and taxonomic level
taxa_counts <- data.frame(
  Database = c("Curada", "Genérica"),
  Genus = c(count_taxa(genus_old)["Genus"], count_taxa(genus_new)["Genus"]),
  Species = c(count_taxa(species_old)["Species"], count_taxa(species_new)["Species"])
)

# Reshape data into long format for plotting
taxa_counts_long <- tidyr::pivot_longer(
  taxa_counts, 
  cols = c("Genus", "Species"), 
  names_to = "TaxonomicLevel", 
  values_to = "Count"
)

# Create bar plot of taxa counts for each database and level
ggplot(taxa_counts_long, aes(x = TaxonomicLevel, y = Count, fill = Database)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  theme_minimal(base_size = 14) +
  labs(
    title = "Number of taxa in each database",
    x = "Taxonomic level",
    y = "Number of taxa",
    fill = "Database"
  ) +
  scale_fill_manual(values = c("Curada" = "#A16928", "Genérica" = "#2887A1")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )
