#!/usr/bin/env Rscript

#' Complete SpatialMetabolics Analysis Pipeline
#'
#' This script demonstrates the full workflow for analyzing metabolic states
#' in Visium HD spatial transcriptomics data comparing wild-type and
#' small vessel disease murine brain samples.
#'
#' @author Your Name
#' @date 2024

# Load required libraries
library(SpatialMetabolics)
library(Seurat)
library(patchwork)
library(viridis)
library(BiocParallel)

# Set up parallel processing
register(MulticoreParam(workers = 4))

# Set parameters
data_dir <- "data/visium_hd/"
output_dir <- "results/"
samples <- list(
  WT1 = "WT_sample1",
  WT2 = "WT_sample2",
  SVD1 = "SVD_sample1",
  SVD2 = "SVD_sample2"
)

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. Load and Create SpatialMetabolic Object
# ============================================================================

message("Loading Visium HD data...")

# Load Seurat objects
seurat_list <- lapply(names(samples), function(s) {
  obj <- Load10X_Spatial(
    data.dir = file.path(data_dir, samples[[s]]),
    filename = "filtered_feature_bc_matrix.h5",
    bin.size = 8  # 8μm resolution for Visium HD
  )

  # Add metadata
  obj$sample <- s
  obj$condition <- ifelse(grepl("WT", s), "WT", "SVD")

  # Basic QC
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)

  return(obj)
})

# Create SpatialMetabolic object
spm <- createSpatialMetabolic(
  seurat_list = seurat_list,
  condition_col = "condition",
  sample_col = "sample",
  min_features = 200,
  min_spots = 10
)

# Save initial object
saveSpatialMetabolic(spm, file.path(output_dir, "spatial_metabolic_raw.rds"))

# ============================================================================
# 2. Quality Control and Preprocessing
# ============================================================================

message("Performing quality control...")

# Calculate QC metrics
spm <- calculateQCMetrics(spm)

# Create QC plots
p_qc <- plotQCSpatial(spm, metrics = c("nCount_RNA", "nFeature_RNA", "percent_mt"))
ggsave(file.path(output_dir, "qc_spatial.pdf"), p_qc, width = 12, height = 8)

# Filter based on QC metrics (if needed)
# spm <- subsetSpatial(spm, cells = which(colData(spm)$percent_mt < 15))

message("Preprocessing spatial data...")

# Normalize and find variable features
spm <- preprocessSpatial(
  spm,
  normalize_method = "logNormalize",
  n_variable_features = 3000,
  n_pcs = 50
)

# Run UMAP
spm <- runUMAP(spm, n_components = 2, use_dims = 30)

# ============================================================================
# 3. Calculate Metabolic Scores
# ============================================================================

message("Calculating metabolic pathway scores...")

# Define all pathways
pathways <- getAllMetabolicPathways("mouse")

# Print pathway coverage
coverage_info <- sapply(pathways, function(genes) {
  present <- sum(genes %in% rownames(spm))
  total <- length(genes)
  percent <- round(present/total * 100, 1)
  return(c(present = present, total = total, percent = percent))
})

write.csv(t(coverage_info),
          file.path(output_dir, "pathway_coverage.csv"),
          row.names = TRUE)

# Calculate scores
spm <- calculateMetabolicScores(
  spm,
  pathways = pathways,
  method = "seurat",
  scale = TRUE,
  nbin = 24,
  BPPARAM = MulticoreParam(workers = 4)
)

# ============================================================================
# 4. Spatial Pattern Detection
# ============================================================================

message("Detecting spatial metabolic gradients...")

spm <- detectMetabolicGradients(
  spm,
  features = NULL,  # Use metabolic scores
  method = "moran",
  n_neighbors = 6,
  permutations = 999,  # Use 999 for publication
  BPPARAM = MulticoreParam(workers = 4)
)

# Extract and save results
gradient_results <- analysisResults(spm)[["spatial_gradients_scores"]]
write.csv(gradient_results,
          file.path(output_dir, "spatial_gradient_results.csv"),
          row.names = FALSE)

# ============================================================================
# 5. Statistical Comparison
# ============================================================================

message("Comparing metabolic states between conditions...")

# Compare pathway scores
pathway_results <- compareMetabolicStates(
  spm,
  condition = "condition",
  ref_group = "WT",
  test_group = "SVD",
  test_type = "wilcox",
  features = NULL,  # Test metabolic scores
  BPPARAM = MulticoreParam(workers = 4)
)

write.csv(pathway_results,
          file.path(output_dir, "pathway_differential_results.csv"),
          row.names = FALSE)

# Compare individual genes in key pathways
key_genes <- unique(c(
  pathways[["OXPHOS"]][1:50],
  pathways[["Glycolysis"]],
  pathways[["ATP_Metabolism"]][1:30]
))

gene_results <- compareMetabolicStates(
  spm,
  condition = "condition",
  ref_group = "WT",
  test_group = "SVD",
  test_type = "wilcox",
  features = intersect(key_genes, rownames(spm)),
  BPPARAM = MulticoreParam(workers = 4)
)

write.csv(gene_results,
          file.path(output_dir, "gene_differential_results.csv"),
          row.names = FALSE)

# Pathway enrichment analysis
enrichment_results <- pathwayEnrichment(
  gene_results,
  pathways = pathways,
  pvalue_cutoff = 0.05,
  log2fc_cutoff = 0.25
)

write.csv(enrichment_results,
          file.path(output_dir, "pathway_enrichment_results.csv"),
          row.names = FALSE)

# ============================================================================
# 6. Find Metabolic Neighborhoods
# ============================================================================

message("Identifying metabolic neighborhoods...")

spm <- findMetabolicNeighborhoods(
  spm,
  k = 5,
  method = "kmeans",
  use_pca = TRUE,
  n_pcs = 10,
  seed = 42
)

# ============================================================================
# 7. Generate Publication Figures
# ============================================================================

message("Creating publication-ready figures...")

# Figure 1: Main results figure
# Panel A: Spatial OXPHOS activity
p1 <- plotSpatialMetabolicScore(
  spm,
  pathway = "OXPHOS",
  samples = c("WT1", "SVD1"),
  ncol = 2,
  point_size = 0.5,
  color_scale = "plasma"
) +
  labs(title = "OXPHOS Activity") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Panel B: Metabolic score comparison
p2 <- plotMetabolicComparison(
  spm,
  pathway = "OXPHOS",
  group_by = "condition",
  plot_type = "violin"
) +
  labs(title = "OXPHOS Comparison") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Panel C: Spatial gradients
p3 <- plotSpatialGradients(
  spm,
  top_n = 10
) +
  labs(title = "Top Spatial Patterns") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Panel D: Pathway heatmap
p4 <- plotMetabolicSummary(
  spm,
  group_by = "condition",
  plot_type = "box",
  show_significance = TRUE
)

# Combine panels
figure1 <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    tag_levels = 'A',
    title = "Spatial Metabolic Dysregulation in Small Vessel Disease",
    subtitle = "Visium HD analysis of murine brain tissue",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# Save figure
ggsave(
  filename = file.path(output_dir, "figure1_main_results.tiff"),
  plot = figure1,
  width = 183,
  height = 200,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = file.path(output_dir, "figure1_main_results.pdf"),
  plot = figure1,
  width = 183,
  height = 200,
  units = "mm"
)

# Figure 2: Comprehensive metabolic analysis
# Spatial grid of all pathways
p_grid <- plotSpatialGrid(
  spm,
  features = names(pathways),
  sample = "WT1",
  ncol = 3,
  sync_scale = FALSE
)

ggsave(
  filename = file.path(output_dir, "figure2_pathway_grid.pdf"),
  plot = p_grid,
  width = 250,
  height = 180,
  units = "mm"
)

# Figure 3: Metabolic transitions
p_scatter <- plotMetabolicScatter(
  spm,
  pathway1 = "OXPHOS",
  pathway2 = "Glycolysis",
  group_by = "condition",
  add_ellipse = TRUE
)

p_volcano <- plotVolcano(
  pathway_results,
  pval_cutoff = 0.05,
  fc_cutoff = 0.25,
  label_top = 10
)

figure3 <- p_scatter | p_volcano

ggsave(
  filename = file.path(output_dir, "figure3_metabolic_shift.pdf"),
  plot = figure3,
  width = 183,
  height = 100,
  units = "mm"
)

# ============================================================================
# 8. Additional Visualizations
# ============================================================================

message("Creating supplementary figures...")

# ATP metabolism spatial plot
p_atp <- plotSpatialMetabolicScore(
  spm,
  pathway = "ATP_Metabolism",
  samples = c("WT1", "WT2", "SVD1", "SVD2"),
  ncol = 2,
  point_size = 0.5,
  color_scale = "inferno"
)

ggsave(
  filename = file.path(output_dir, "supp_figure_atp_spatial.pdf"),
  plot = p_atp,
  width = 183,
  height = 183,
  units = "mm"
)

# Metabolic neighborhoods
p_neighborhoods <- plotSpatialFeatures(
  spm,
  features = "metabolic_neighborhood",
  ncol = 2,
  point_size = 0.5
)

ggsave(
  filename = file.path(output_dir, "supp_figure_neighborhoods.pdf"),
  plot = p_neighborhoods,
  width = 183,
  height = 120,
  units = "mm"
)

# Sample correlation
p_cor <- plotSampleCorrelation(spm, use_pathways = TRUE)

# Network analysis (if enough pathways)
if (nrow(metabolicScores(spm)) >= 5) {
  p_network <- plotMetabolicNetwork(spm, min_correlation = 0.5)

  ggsave(
    filename = file.path(output_dir, "supp_figure_network.pdf"),
    plot = p_network,
    width = 150,
    height = 150,
    units = "mm"
  )
}

# ============================================================================
# 9. Generate Analysis Report
# ============================================================================

message("Generating analysis report...")

report <- createMetabolicReport(
  spm,
  save_to = file.path(output_dir, "metabolic_analysis_report.pdf")
)

# ============================================================================
# 10. Save Final Results
# ============================================================================

message("Saving final analysis results...")

# Save processed SpatialMetabolic object
saveSpatialMetabolic(spm, file.path(output_dir, "spatial_metabolic_processed.rds"))

# Create results summary
summary_stats <- list(
  n_spots_total = ncol(spm),
  n_spots_per_sample = table(colData(spm)$sample_id),
  n_features = nrow(spm),
  n_variable_features = sum(rowData(spm)$variable_feature),
  pathway_scores_summary = summary(metabolicScores(spm)),
  significant_pathways = pathway_results[pathway_results$adj_pvalue < 0.05, ],
  top_differential_genes = head(gene_results[order(gene_results$adj_pvalue), ], 20),
  spatial_patterns = head(gradient_results[order(gradient_results$adj_pvalue), ], 10),
  enriched_pathways = enrichment_results[enrichment_results$adj_pvalue < 0.1, ]
)

saveRDS(
  summary_stats,
  file = file.path(output_dir, "analysis_summary.rds")
)

# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "session_info.txt")
)

message("\nAnalysis complete! Results saved to: ", output_dir)
message("\nKey findings:")
message("- ", nrow(summary_stats$significant_pathways), " significantly altered pathways")
message("- ", sum(gene_results$adj_pvalue < 0.05), " differentially expressed metabolic genes")
message("- ", sum(gradient_results$adj_pvalue < 0.05), " pathways with spatial patterns")

# Print top results
message("\nTop altered pathways:")
print(head(pathway_results[order(pathway_results$adj_pvalue), c("feature", "log2FC", "adj_pvalue")], 5))
