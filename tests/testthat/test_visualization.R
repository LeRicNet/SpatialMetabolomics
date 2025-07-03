# Test visualization functions

test_that("Spatial metabolic score plotting works", {
  # Create test data with scores
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  # Use genes that actually exist in the test data
  pathways <- list(TestPathway = c("Gene1", "Gene2", "Gene3"))
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  # Basic plot
  p <- plotSpatialMetabolicScore(spm, pathway = "TestPathway")
  expect_ggplot(p)

  # Test different parameters
  p2 <- plotSpatialMetabolicScore(
    spm,
    pathway = "TestPathway",
    samples = "Sample1",
    color_scale = "plasma",
    point_size = 2
  )
  expect_ggplot(p2)

  # Check error handling
  expect_error(
    plotSpatialMetabolicScore(spm, pathway = "NonExistent"),
    "not found in metabolic scores"
  )
})

test_that("Metabolic comparison plots work", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  pathways <- list(TestPathway = c("Gene1", "Gene2", "Gene3"))
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  # Violin plot
  p_violin <- plotMetabolicComparison(
    spm,
    pathway = "TestPathway",
    plot_type = "violin",
    add_stats = FALSE
  )
  expect_ggplot(p_violin)

  # Box plot
  p_box <- plotMetabolicComparison(
    spm,
    pathway = "TestPathway",
    plot_type = "box"
  )
  expect_ggplot(p_box)

  # Density plot
  p_density <- plotMetabolicComparison(
    spm,
    pathway = "TestPathway",
    plot_type = "density"
  )
  expect_ggplot(p_density)
})

test_that("Spatial gradient plotting works", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  pathways <- list(TestPathway = c("Gene1", "Gene2", "Gene3"))
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  # Detect gradients first
  spm <- detectMetabolicGradients(spm, method = "moran", permutations = 0)

  # Plot gradients
  p <- plotSpatialGradients(spm, top_n = 5)
  expect_ggplot(p)

  # Check data structure
  expect_true("feature" %in% names(p$data))
})

test_that("Volcano plot works", {
  # Create DE results
  de_results <- data.frame(
    feature = paste0("Gene", 1:100),
    log2FC = rnorm(100, sd = 1),
    pvalue = runif(100),
    adj_pvalue = p.adjust(runif(100)),
    stringsAsFactors = FALSE
  )

  # Basic volcano
  p <- plotVolcano(de_results, label_top = 0)
  expect_ggplot(p)

  # With labels (need ggrepel)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_labeled <- plotVolcano(de_results, label_top = 5)
    expect_ggplot(p_labeled)
  }

  # Custom thresholds
  p_custom <- plotVolcano(
    de_results,
    pval_cutoff = 0.01,
    fc_cutoff = 1,
    label_top = 0
  )
  expect_ggplot(p_custom)
})

test_that("Heatmap plotting works", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  pathways <- list(
    Pathway1 =  c("Gene1", "Gene2", "Gene3"),
    Pathway2 = c("Gene4", "Gene5", "Gene6")
  )
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  # Test with metabolic scores
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    p <- plotMetabolicHeatmap(spm, scale_rows = TRUE)
    expect_true(!is.null(p))
  } else {
    # ggplot2 fallback
    p <- plotMetabolicHeatmap(spm, scale_rows = TRUE)
    expect_ggplot(p)
  }

  # Test with specific features
  p2 <- plotMetabolicHeatmap(
    spm,
    features = c("Gene1", "Gene3"),
    show_row_names = TRUE,
    scale_rows = FALSE
  )
  expect_true(!is.null(p2))
})

test_that("Spatial feature plotting works", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  pathways <- list(TestPathway = c("Gene1", "Gene2", "Gene3"))
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  # Plot genes
  p_genes <- plotSpatialFeatures(
    spm,
    features = c("Gene1", "Gene2"),
    point_size = 1
  )
  expect_ggplot(p_genes)

  # Plot scores
  p_scores <- plotSpatialFeatures(
    spm,
    features = "TestPathway",
    sync_scales = TRUE
  )
  expect_ggplot(p_scores)

  # Mixed features
  p_mixed <- plotSpatialFeatures(
    spm,
    features = c("Gene1", "TestPathway"),
    ncol = 2
  )
  expect_ggplot(p_mixed)
})

test_that("Summary plots work", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  pathways <- list(
    Pathway1 = c("Gene1", "Gene2", "Gene3"),
    Pathway2 = c("Gene4", "Gene5", "Gene6")
  )
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  # Summary plot
  p <- plotMetabolicSummary(
    spm,
    group_by = "condition",
    plot_type = "box",
    show_significance = FALSE
  )
  expect_ggplot(p)
})

test_that("Color scales work correctly", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  pathways <- list(TestPathway = c("Gene1", "Gene2", "Gene3"))
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  # Test different color scales
  color_scales <- c("viridis", "plasma", "inferno", "magma", "cividis")

  for (cs in color_scales) {
    p <- plotSpatialMetabolicScore(
      spm,
      pathway = "TestPathway",
      color_scale = cs
    )
    expect_ggplot(p)
  }
})

test_that("Plot error handling works", {
  spm <- create_test_data()

  # No scores calculated
  expect_error(
    plotSpatialMetabolicScore(spm, "Test"),
    "No metabolic scores calculated"
  )

  expect_error(
    plotMetabolicComparison(spm, "Test"),
    "not found"
  )

  # Invalid parameters
  spm <- normalizeSpatial(spm, verbose = FALSE)
  pathways <- list(TestPathway = c("Gene1", "Gene2", "Gene3"))
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  expect_error(
    plotMetabolicComparison(spm, "TestPathway", group_by = "invalid_column"),
    "not found in colData"
  )
})

test_that("Additional visualization functions work", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  spm <- calculateQCMetrics(spm)
  pathways <- list(
    Pathway1 = c("Gene1", "Gene2", "Gene3"),
    Pathway2 = c("Gene4", "Gene5", "Gene6")
  )
  spm <- calculateMetabolicScores(spm, pathways = pathways)  # REMOVED verbose = FALSE

  # QC spatial plot
  p_qc <- plotQCSpatial(spm, metrics = c("nCount_RNA", "nFeature_RNA"))
  expect_ggplot(p_qc)

  # Scatter plot
  p_scatter <- plotMetabolicScatter(
    spm,
    pathway1 = "Pathway1",
    pathway2 = "Pathway2",
    add_ellipse = FALSE
  )
  expect_ggplot(p_scatter)

  # Sample correlation (returns pheatmap or ggplot)
  p_cor <- plotSampleCorrelation(spm, use_pathways = TRUE)
  expect_true(!is.null(p_cor))
})
