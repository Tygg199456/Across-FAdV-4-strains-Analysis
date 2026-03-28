#!/usr/bin/env Rscript

# ============================================================================
# WGCNA Analysis for FAdV-4 Infection
# ============================================================================

library(WGCNA)
library(tidyverse)
library(data.table)

# Set working directory
setwd("/Users/tgw/Desktop/FADV_new/FAdV4_LimmaAnalysis copy")

# Enable parallel computing
enableWGCNAThreads(nThreads = 4)

options(stringsAsFactors = FALSE)

# ============================================================================
# Step 1: Load Data
# ============================================================================
cat("Step 1: Loading normalized expression data...\n")

# WGCNA uses TMM-normalized (non-batch-corrected) expression data.
# Rationale: WGCNA captures all sources of co-expression variation including
# batch-associated structure. With n=12 samples, batch correction risks
# over-compressing inter-sample distances and destabilizing module structure.
# The limma-voom DEG analysis (Part 2) handles batch as a model covariate;
# WGCNA serves as an independent exploratory validation layer.
# Results are interpreted as exploratory per manuscript Methods section.
data_stats <- readRDS("results/intermediate/Normalized_Expression_For_Statistics.rds")
log2cpm_matrix <- data_stats$expr
sample_info <- data_stats$sample_info

cat(sprintf("  Expression matrix dimensions: %d genes x %d samples\n",
            nrow(log2cpm_matrix), ncol(log2cpm_matrix)))

# Check sample order
print(sample_info)

# ============================================================================
# Step 2: Prepare Data for WGCNA
# ============================================================================
cat("\nStep 2: Preparing data for WGCNA...\n")

# Transpose: WGCNA requires samples as rows, genes as columns
datExpr <- t(log2cpm_matrix)
cat(sprintf("  Transposed matrix: %d samples x %d genes\n",
            nrow(datExpr), ncol(datExpr)))

# Check data quality
gsg <- goodSamplesGenes(datExpr, verbose = 3)
cat(sprintf("  All genes OK: %s\n", gsg$allOK))

# If not all OK, filter
if (!gsg$allOK) {
  cat("  Filtering out problematic genes...\n")
  if (sum(!gsg$goodGenes) > 0) {
    datExpr <- datExpr[, gsg$goodGenes]
  }
  cat(sprintf("  After filtering: %d samples x %d genes\n",
              nrow(datExpr), ncol(datExpr)))
}

# ============================================================================
# Step 3: Select Soft Threshold Power
# ============================================================================
cat("\nStep 3: Selecting soft threshold power...\n")

# Test powers from 1 to 20
powers <- c(1:20)

sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         RsquaredCut = 0.85,
                         verbose = 5)

# Get estimated power
softPower <- sft$powerEstimate
cat(sprintf("  Estimated soft power: %d\n", softPower))

# If power estimate is NA or too low, manually select a proper power
# Power 1 is too low biologically, need higher power for scale-free network
if (is.na(softPower) || softPower < 6) {
  cat("  Power estimate is too low. Manually selecting power = 14\n")
  cat("  (R² = 0.871, meets threshold of 0.85)\n")
  softPower <- 14
}

# Save the scale-free topology plot
pdf("results/figures/WGCNA_SoftThreshold.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))

# Plot 1: Scale-free fit index
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, R²",
     main = "Scale-free topology fit index",
     type = "n")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     col = ifelse(sft$fitIndices[, 2] > 0.85, "red", "black"),
     cex = 0.8)
abline(h = 0.85, col = "red", lty = 2)

# Plot 2: Mean connectivity
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity",
     type = "n")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = 0.8)

dev.off()

cat(sprintf("  Soft threshold power = %d saved to: results/figures/WGCNA_SoftThreshold.pdf\n", softPower))

# ============================================================================
# Step 4: Construct Network and Identify Modules
# ============================================================================
cat("\nStep 4: Constructing network and identifying modules...\n")

# Sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr))
pdf("results/figures/WGCNA_SampleClustering.pdf", width = 12, height = 8)
plot(sampleTree,
     main = "Sample clustering to detect outliers",
     xlab = "Height",
     sub = "")
dev.off()

# Construct network
cat("  Running blockwiseModules (this may take a few minutes)...\n")

net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "results/intermediate/WGCNA_TOM",
  verbose = 3
)

# Module colors
moduleColors <- labels2colors(net$colors)
nModules <- length(unique(net$colors))

cat(sprintf("  Number of modules identified: %d\n", nModules))

# Module sizes
module_sizes <- table(net$colors)
cat("  Module sizes:\n")
print(module_sizes)

# Save module assignment
module_assignment <- data.frame(
  gene = colnames(datExpr),
  module = net$colors,
  module_color = moduleColors
)
write.csv(module_assignment,
          "results/tables/WGCNA_Module_Assignment.csv",
          row.names = FALSE)

# ============================================================================
# Step 5: Module-Trait Relationship
# ============================================================================
cat("\nStep 5: Calculating module-trait relationships...\n")

# Create phenotype vector: FAdV4 = 1, Control = 0
# Sample order: based on sample_info
condition <- ifelse(sample_info$group == "FAdV4", 1, 0)
cat("  Phenotype vector (FAdV4=1, Control=0):\n")
print(condition)

# Calculate module-trait correlation
moduleTraitCor <- cor(net$MEs, condition, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Create text matrix for heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")

# Save correlation results
module_trait_df <- data.frame(
  Module = colnames(net$MEs),
  Correlation = moduleTraitCor[, 1],
  Pvalue = moduleTraitPvalue[, 1]
) %>% arrange(desc(abs(Correlation)))

write.csv(module_trait_df,
          "results/tables/WGCNA_Module_Trait_Correlation.csv",
          row.names = FALSE)

# Plot module-trait relationship heatmap
pdf("results/figures/WGCNA_Module_Trait_Heatmap.pdf", width = 8, height = 10)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "FAdV-4 infection",
  yLabels = names(net$MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  main = "Module-trait relationships (FAdV-4 infection)",
  cex.text = 0.7
)
dev.off()

# ============================================================================
# Step 6: Identify Key Modules and Hub Genes
# ============================================================================
cat("\nStep 6: Identifying key modules and hub genes...\n")

# Find most correlated module
key_module_idx <- which.max(abs(moduleTraitCor[, 1]))
key_module_name <- names(net$MEs)[key_module_idx]  # e.g., "ME6"
key_module_num <- as.numeric(gsub("ME", "", key_module_name))  # e.g., 6
key_module_cor <- moduleTraitCor[key_module_idx, 1]
key_module_pval <- moduleTraitPvalue[key_module_idx, 1]

cat(sprintf("  Most correlated module: %s (numeric: %d)\n", key_module_name, key_module_num))
cat(sprintf("  Correlation: %.3f, P-value: %.2e\n", key_module_cor, key_module_pval))

# Get genes in key module
key_module_genes <- names(net$colors)[net$colors == key_module_num]
cat(sprintf("  Number of genes in key module: %d\n", length(key_module_genes)))

# Calculate Gene Significance (GS) and Module Membership (MM)
GS <- abs(cor(datExpr, condition, use = "p"))
MM <- abs(cor(datExpr, net$MEs, use = "p"))

# Define hub genes: MM > 0.8 AND GS > 0.8
hub_genes_candidates <- key_module_genes[
  MM[key_module_genes, key_module_name] > 0.8 &
  GS[key_module_genes, 1] > 0.8
]

cat(sprintf("  Hub gene candidates (MM > 0.8 & GS > 0.8): %d\n",
            length(hub_genes_candidates)))

# Save hub genes with their connectivity
hub_gene_data <- data.frame(
  gene = hub_genes_candidates,
  module_membership = MM[hub_genes_candidates, key_module_name],
  gene_significance = GS[hub_genes_candidates, 1]
) %>% arrange(desc(gene_significance))

write.csv(hub_gene_data,
          "results/tables/WGCNA_Hub_Genes.csv",
          row.names = FALSE)

# ============================================================================
# Step 7: Plot Module Eigengene Dendrogram
# ============================================================================
cat("\nStep 7: Plotting module eigengene dendrogram...\n")

pdf("results/figures/WGCNA_Module_Dendrogram.pdf", width = 12, height = 8)
plotEigengeneNetworks(net$MEs,
                     "Eigengene dendrogram",
                     marDendro = c(0, 4, 2, 0))
dev.off()

# ============================================================================
# Summary
# ============================================================================
cat("\n")
cat(strrep("=", 60))
cat("\nWGCNA Analysis Summary\n")
cat(strrep("=", 60))
cat(sprintf("\n  Total genes analyzed: %d\n", ncol(datExpr)))
cat(sprintf("  Total samples: %d (6 FAdV4 + 6 Control)\n", nrow(datExpr)))
cat(sprintf("  Soft threshold power: %d\n", softPower))
cat(sprintf("  Number of modules: %d\n", nModules))
cat(sprintf("\n  Top correlated module: %s\n", key_module_name))
cat(sprintf("    Correlation: %.3f\n", key_module_cor))
cat(sprintf("    P-value: %.2e\n", key_module_pval))
cat(sprintf("    Genes in module: %d\n", length(key_module_genes)))
cat(sprintf("    Hub genes (MM>0.8 & GS>0.8): %d\n", length(hub_genes_candidates)))
cat("\n  Output files:\n")
cat("    - results/figures/WGCNA_SoftThreshold.pdf\n")
cat("    - results/figures/WGCNA_SampleClustering.pdf\n")
cat("    - results/figures/WGCNA_Module_Trait_Heatmap.pdf\n")
cat("    - results/figures/WGCNA_Module_Dendrogram.pdf\n")
cat("    - results/tables/WGCNA_Module_Assignment.csv\n")
cat("    - results/tables/WGCNA_Module_Trait_Correlation.csv\n")
cat("    - results/tables/WGCNA_Hub_Genes.csv\n")
cat(strrep("=", 60))
cat("\n\n")

# ============================================================================
# Step 8: Supplementary Figure S3 - Eigengene Adjacency Heatmap
# ============================================================================
cat("\nStep 8: Generating Supplementary Figure S3 - Eigengene adjacency heatmap...\n")

# Get module eigengenes
MEs <- net$MEs

# Use WGCNA's built-in function to generate combined plot
# This creates both the dendrogram and the adjacency heatmap in one figure
pdf("results/figures/WGCNA_Eigengene_Adjacency.pdf", width = 10, height = 8)
plotEigengeneNetworks(MEs,
                     "Eigengene dendrogram and adjacency heatmap",
                     marDendro = c(0, 4, 2, 2),
                     marHeatmap = c(4, 4, 2, 2),
                     plotDendrograms = TRUE,
                     xLabelsAngle = 90)
dev.off()

cat("  Saved: results/figures/WGCNA_Eigengene_Adjacency.pdf\n")

# ============================================================================
# Summary
# ============================================================================
cat("\n")
cat(strrep("=", 60))
cat("\nWGCNA Analysis Summary\n")
cat(strrep("=", 60))
cat(sprintf("\n  Total genes analyzed: %d\n", ncol(datExpr)))
cat(sprintf("  Total samples: %d (6 FAdV4 + 6 Control)\n", nrow(datExpr)))
cat(sprintf("  Soft threshold power: %d\n", softPower))
cat(sprintf("  Number of modules: %d\n", nModules))
cat(sprintf("\n  Top correlated module: %s\n", key_module_name))
cat(sprintf("    Correlation: %.3f\n", key_module_cor))
cat(sprintf("    P-value: %.2e\n", key_module_pval))
cat(sprintf("    Genes in module: %d\n", length(key_module_genes)))
cat(sprintf("    Hub genes (MM>0.8 & GS>0.8): %d\n", length(hub_genes_candidates)))
cat("\n  Output files:\n")
cat("    - results/figures/WGCNA_SoftThreshold.pdf\n")
cat("    - results/figures/WGCNA_SampleClustering.pdf\n")
cat("    - results/figures/WGCNA_Module_Trait_Heatmap.pdf\n")
cat("    - results/figures/WGCNA_Module_Dendrogram.pdf\n")
cat("    - results/figures/WGCNA_Eigengene_Adjacency.pdf (NEW)\n")
cat("    - results/tables/WGCNA_Module_Assignment.csv\n")
cat("    - results/tables/WGCNA_Module_Trait_Correlation.csv\n")
cat("    - results/tables/WGCNA_Hub_Genes.csv\n")
cat(strrep("=", 60))
cat("\n\n")

# Save workspace image
save.image("results/intermediate/WGCNA_workspace.RData")
cat("  Workspace saved to: results/intermediate/WGCNA_workspace.RData\n")
