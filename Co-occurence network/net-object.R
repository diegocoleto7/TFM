library(phyloseq)
library(microbiome)
library(genefilter)

# Read RDS (phyloseq) objects
AI_f <- readRDS("C:/Users/Diego/Desktop/TFM/Redes-Bacteria-Hongo/Red1/phyloseq_WGS_AI.rds")
Celiac_f <- readRDS("C:/Users/Diego/Desktop/TFM/Redes-Bacteria-Hongo/Red1/phyloseq_WGS_celiacos.rds")
AI_b <- readRDS("C:/Users/Diego/Desktop/TFM/Redes-Bacteria-Hongo/physeq_mpa3_v1_bac.rds")
Celiac_b <- readRDS("C:/Users/Diego/Desktop/TFM/Redes-Bacteria-Hongo/physeq_celiac_bac.rds")

# Taxonomic agglomeration for Celiac bacteria at Species level
Celiac_b <- tax_glom(Celiac_b, taxrank = "Species")

# Filtering based on relative abundance
flist_celiac <- filterfun(kOverA(1, 0.05))  # Keep taxa > 0.05% in ≥ 1 sample (Celiac_b)
Celiac_b.fil <- filter_taxa(Celiac_b, flist_celiac)
Celiac_b <- prune_taxa(Celiac_b.fil, Celiac_b)

flist_AI <- filterfun(kOverA(k = 2, A = 0.005))  # Keep taxa > 0.005% in ≥ 2 samples (AI_b)
AI_b.fill <- filter_taxa(AI_b, flist_AI)
AI_b <- prune_taxa(AI_b.fill, AI_b)

# Transform counts to parts-per-million (PPM)
AI_b <- transform_sample_counts(AI_b, function(x) 1E6 * x)
Celiac_b <- transform_sample_counts(Celiac_b, function(x) 1E6 * x)

# Remove phylogenetic trees from phyloseq objects
AI_f@phy_tree <- NULL 
AI_b@phy_tree <- NULL 
Celiac_f@phy_tree <- NULL 
Celiac_b@phy_tree <- NULL

# Filter samples by matching sample names (Celiac)
otu_cols_Celiac_f <- colnames(otu_table(Celiac_f))
otu_cols_Celiac_b <- colnames(otu_table(Celiac_b))
only_in_f_celiac <- setdiff(otu_cols_Celiac_f, otu_cols_Celiac_b)
only_in_b_celiac <- setdiff(otu_cols_Celiac_b, otu_cols_Celiac_f)
in_both_celiac <- intersect(otu_cols_Celiac_f, otu_cols_Celiac_b)
cat("Samples only in Celiac_f:\n", only_in_f_celiac, "\n")
cat("Samples only in Celiac_b:\n", only_in_b_celiac, "\n")
cat("Samples in both Celiac_f & Celiac_b:\n", in_both_celiac, "\n")

# Filter samples by matching sample names (AI)
otu_cols_AI_f <- colnames(otu_table(AI_f))
otu_cols_AI_b <- colnames(otu_table(AI_b))
only_in_f <- setdiff(otu_cols_AI_f, otu_cols_AI_b)
only_in_b <- setdiff(otu_cols_AI_b, otu_cols_AI_f)
in_both <- intersect(otu_cols_AI_f, otu_cols_AI_b)
cat("Samples only in AI_f:\n", only_in_f, "\n")
cat("Samples only in AI_b:\n", only_in_b, "\n")
cat("Samples in both AI_f & AI_b:\n", in_both, "\n")

# Prune phyloseq objects to keep only matching samples for AI
AI_f <- prune_samples(sample_names(AI_f) %in% in_both, AI_f)
AI_b <- prune_samples(sample_names(AI_b) %in% in_both, AI_b)

# Rename taxonomic labels
tax_df_celiac <- as.data.frame(tax_table(Celiac_f))
new_taxa_names <- paste(tax_df_celiac$Genus, tax_df_celiac$Species, sep = "_")
taxa_names(Celiac_f) <- new_taxa_names

tax_df_AI <- as.data.frame(tax_table(AI_f))
new_taxa_names <- paste(tax_df_AI$Genus, tax_df_AI$Species, sep = "_")
taxa_names(AI_f) <- new_taxa_names

taxa_names(AI_b) <- tax_table(AI_b)[, "Species"]

# Merge fungi and bacteria phyloseq objects to create integrated objects
AI_net <- merge_phyloseq(AI_b, AI_f)
Celiac_net <- merge_phyloseq(Celiac_b, Celiac_f)

# Save final phyloseq objects for network analyses
saveRDS(AI_net, file = "C:/Users/Diego/Desktop/TFM/Redes-Bacteria-Hongo/Red1/AI-net-diego.rds")
saveRDS(Celiac_net, file = "C:/Users/Diego/Desktop/TFM/Redes-Bacteria-Hongo/Red1/Celiac-net-diego.rds")
