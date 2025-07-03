# Unit tests for SpatialMetabolic package

library(testthat)
library(SpatialMetabolics)
library(SpatialExperiment)
library(SingleCellExperiment)

# Helper function to create test data
create_test_data <- function(n_genes = 100, n_spots = 50, n_samples = 2) {
  # Create count matrices
  counts_list <- lapply(seq_len(n_samples), function(i) {
    matrix(
      rpois(n_genes * n_spots, lambda = 5),
      nrow = n_genes,
      ncol = n_spots,
      dimnames = list(
        paste0("Gene", seq_len(n_genes)),
        paste0("Spot", seq_len(n_spots), "_S", i)
      )
    )
  })

  # Combine counts
  counts <- do.call(cbind, counts_list)

  # Create metadata
  col_data <- data.frame(
    sample_id = rep(paste0("Sample", seq_len(n_samples)), each = n_spots),
    condition = rep(c("WT", "Disease"), each = n_spots)[seq_len(n_spots * n_samples)],
    row.names = colnames(counts)
  )

  # Create spatial coordinates
  coords <- do.call(rbind, lapply(seq_len(n_samples), function(i) {
    data.frame(
      x = rep(1:10, 5),
      y = rep(1:5, each = 10)
    )
  }))

  # Create SpatialExperiment
  spe <- SpatialExperiment(
    assays = list(counts = counts),
    colData = col_data,
    spatialCoords = as.matrix(coords)
  )

  # Convert to SpatialMetabolic
  spm <- as(spe, "SpatialMetabolic")
  comparisonGroups(spm) <- factor(col_data$condition)

  return(spm)
}

# Test class creation and validity
test_that("SpatialMetabolic class creation works", {
  spm <- create_test_data()

  expect_s4_class(spm, "SpatialMetabolic")
  expect_true(inherits(spm, "SpatialExperiment"))
  expect_equal(ncol(spm), 100)
  expect_equal(nrow(spm), 100)
})

# Test accessor methods
test_that("Accessor methods work correctly", {
  spm <- create_test_data()

  # Test metabolic pathways accessor
  test_pathways <- list(
    pathway1 = c("Gene1", "Gene2", "Gene3"),
    pathway2 = c("Gene4", "Gene5", "Gene6")
  )

  metabolicPathways(spm) <- test_pathways
  expect_equal(metabolicPathways(spm), test_pathways)

  # Test metabolic scores accessor
  test_scores <- matrix(rnorm(2 * ncol(spm)), nrow = 2, ncol = ncol(spm))
  rownames(test_scores) <- names(test_pathways)

  metabolicScores(spm) <- test_scores
  expect_equal(metabolicScores(spm), test_scores)

  # Test comparison groups
  expect_equal(levels(comparisonGroups(spm)), c("WT", "Disease"))
})

# Test normalization
test_that("Normalization works correctly", {
  spm <- create_test_data()

  # Test log normalization
  spm_norm <- normalizeSpatial(spm, method = "logNormalize", verbose = FALSE)

  expect_true(all(logcounts(spm_norm) >= 0))
  expect_true(max(logcounts(spm_norm)) > 0)

  # Check that normalized values are different from counts
  expect_false(all(counts(spm_norm) == logcounts(spm_norm)))
})

# Test variable feature selection
test_that("Variable feature selection works", {
  spm <- create_test_data(n_genes = 200)
  spm <- normalizeSpatial(spm, verbose = FALSE)

  # Select variable features
  spm <- findVariableFeatures(spm, n_features = 50, verbose = FALSE)

  expect_true("variable_feature" %in% colnames(rowData(spm)))
  expect_equal(sum(rowData(spm)$variable_feature), 50)
})

# Test PCA
test_that("PCA calculation works", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)
  spm <- findVariableFeatures(spm, n_features = 50, verbose = FALSE)

  # Run PCA
  spm <- runPCA(spm, n_components = 10, verbose = FALSE)

  expect_true("PCA" %in% reducedDimNames(spm))
  expect_equal(ncol(reducedDim(spm, "PCA")), 10)
  expect_equal(nrow(reducedDim(spm, "PCA")), ncol(spm))
})

# Test metabolic pathway functions
test_that("Metabolic pathway retrieval works", {
  # Test OXPHOS genes
  oxphos_mouse <- getOXPHOSGenes("mouse")
  expect_true(length(oxphos_mouse) > 100)
  expect_true(all(grepl("^[A-Z][a-z]", oxphos_mouse)))  # Mouse gene format

  # Test glycolysis genes
  glycolysis_mouse <- getGlycolysisGenes("mouse")
  expect_true(length(glycolysis_mouse) > 15)
  expect_true("Hk1" %in% glycolysis_mouse)

  # Test human genes
  oxphos_human <- getOXPHOSGenes("human")
  expect_true(all(grepl("^[A-Z]+[0-9]*", oxphos_human)))  # Human gene format

  # Test pathway list
  all_pathways <- getAllMetabolicPathways("mouse")
  expect_equal(length(all_pathways), 6)
  expect_true("OXPHOS" %in% names(all_pathways))
})

# Test metabolic score calculation
test_that("Metabolic score calculation works", {
  # Create data with known pathway genes
  genes <- c(getOXPHOSGenes("mouse")[1:20],
             getGlycolysisGenes("mouse")[1:10],
             paste0("RandomGene", 1:70))

  spm <- create_test_data(n_genes = 100)
  rownames(spm) <- genes[1:100]
  spm <- normalizeSpatial(spm, verbose = FALSE)

  # Calculate scores
  spm <- calculateMetabolicScores(
    spm,
    pathways = c("oxphos", "glycolysis"),
    method = "seurat",
    scale = FALSE,
    BPPARAM = SerialParam()
  )

  expect_equal(nrow(metabolicScores(spm)), 2)
  expect_equal(ncol(metabolicScores(spm)), ncol(spm))
  expect_true(all(c("OXPHOS", "Glycolysis") %in% rownames(metabolicScores(spm))))

  # Check that scores are added to colData
  expect_true("score_OXPHOS" %in% colnames(colData(spm)))
  expect_true("score_Glycolysis" %in% colnames(colData(spm)))
})

# Test spatial gradient detection
test_that("Spatial gradient detection works", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)

  # Create artificial spatial pattern
  coords <- spatialCoords(spm)
  spatial_feature <- sin(coords[, 1] / 5) + cos(coords[, 2] / 5)
  logcounts(spm)["Gene1", ] <- spatial_feature + rnorm(ncol(spm), sd = 0.1)

  # Calculate metabolic scores first
  pathways <- list(TestPathway = c("Gene1", "Gene2", "Gene3"))
  spm <- calculateMetabolicScores(spm, pathways = pathways, verbose = FALSE)

  # Detect gradients
  spm <- detectMetabolicGradients(
    spm,
    method = "moran",
    n_neighbors = 4,
    permutations = 0,  # Skip permutation for speed
    BPPARAM = SerialParam()
  )

  expect_true("spatial_gradients_scores" %in% names(analysisResults(spm)))

  results <- analysisResults(spm)[["spatial_gradients_scores"]]
  expect_s3_class(results, "data.frame")
  expect_true("morans_i" %in% colnames(results))
  expect_true("pvalue" %in% colnames(results))
})

# Test differential expression
test_that("Metabolic state comparison works", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)

  # Add some differential expression
  logcounts(spm)[1:10, comparisonGroups(spm) == "Disease"] <-
    logcounts(spm)[1:10, comparisonGroups(spm) == "Disease"] + 2

  # Compare states
  results <- compareMetabolicStates(
    spm,
    condition = "condition",
    ref_group = "WT",
    test_group = "Disease",
    test_type = "wilcox",
    features = rownames(spm)[1:20],
    verbose = FALSE,
    BPPARAM = SerialParam()
  )

  expect_s3_class(results, "data.frame")
  expect_true(all(c("feature", "log2FC", "pvalue", "adj_pvalue") %in% colnames(results)))
  expect_equal(nrow(results), 20)

  # Check that some genes are significant
  expect_true(any(results$pvalue < 0.05))
})

# Test visualization functions
test_that("Basic visualization functions work", {
  spm <- create_test_data()
  spm <- normalizeSpatial(spm, verbose = FALSE)

  # Calculate scores for plotting
  pathways <- list(TestPathway = paste0("Gene", 1:10))
  spm <- calculateMetabolicScores(spm, pathways = pathways, verbose = FALSE)

  # Test spatial plot
  p1 <- plotSpatialMetabolicScore(spm, pathway = "TestPathway")
  expect_s3_class(p1, "ggplot")

  # Test comparison plot
  p2 <- plotMetabolicComparison(spm, pathway = "TestPathway", add_stats = FALSE)
  expect_s3_class(p2, "ggplot")

  # Test volcano plot
  de_results <- data.frame(
    feature = paste0("Gene", 1:100),
    log2FC = rnorm(100),
    adj_pvalue = runif(100),
    stringsAsFactors = FALSE
  )
  p3 <- plotVolcano(de_results, label_top = 0)
  expect_s3_class(p3, "ggplot")
})

# Test utility functions
test_that("Utility functions work correctly", {
  spm <- create_test_data()

  # Test QC metrics
  spm <- calculateQCMetrics(spm)
  expect_true(all(c("nCount_RNA", "nFeature_RNA", "percent_mt") %in% colnames(colData(spm))))

  # Test subsetting
  spm_sub <- subsetSpatial(spm, cells = 1:20, features = 1:50)
  expect_equal(ncol(spm_sub), 20)
  expect_equal(nrow(spm_sub), 50)

  # Test saving and loading
  temp_file <- tempfile(fileext = ".rds")
  saveSpatialMetabolic(spm, temp_file)
  expect_true(file.exists(temp_file))

  spm_loaded <- loadSpatialMetabolic(temp_file)
  expect_s4_class(spm_loaded, "SpatialMetabolic")

  unlink(temp_file)
})

# Test error handling
test_that("Functions handle errors appropriately", {
  spm <- create_test_data()

  # Test error when no scores calculated
  expect_error(
    plotSpatialMetabolicScore(spm, pathway = "OXPHOS"),
    "No metabolic scores calculated"
  )

  # Test error with invalid pathway
  spm <- normalizeSpatial(spm, verbose = FALSE)
  spm <- calculateMetabolicScores(spm, pathways = list(Test = "Gene1"), verbose = FALSE)
  expect_error(
    plotSpatialMetabolicScore(spm, pathway = "InvalidPathway"),
    "not found in metabolic scores"
  )

  # Test error with invalid species
  expect_error(
    getOXPHOSGenes("invalid_species"),
    "Species must be"
  )
})

# Test edge cases
test_that("Functions handle edge cases", {
  # Empty object
  spm_empty <- create_test_data(n_genes = 0, n_spots = 0)
  expect_equal(nrow(spm_empty), 0)
  expect_equal(ncol(spm_empty), 0)

  # Single cell/gene
  spm_single <- create_test_data(n_genes = 1, n_spots = 1)
  spm_single <- normalizeSpatial(spm_single, verbose = FALSE)
  expect_equal(dim(logcounts(spm_single)), c(1, 1))

  # All zero counts
  spm_zero <- create_test_data()
  counts(spm_zero)[] <- 0
  spm_zero <- normalizeSpatial(spm_zero, verbose = FALSE)
  expect_true(all(logcounts(spm_zero) == 0))
})

# Test integration between functions
test_that("Full workflow integration works", {
  # Create data
  spm <- create_test_data(n_genes = 200, n_spots = 100)

  # Add real gene names for some genes
  gene_names <- rownames(spm)
  gene_names[1:30] <- c(getOXPHOSGenes("mouse")[1:20], getGlycolysisGenes("mouse")[1:10])
  rownames(spm) <- gene_names

  # Full preprocessing
  spm <- normalizeSpatial(spm, verbose = FALSE)
  spm <- findVariableFeatures(spm, n_features = 100, verbose = FALSE)
  spm <- runPCA(spm, n_components = 20, verbose = FALSE)

  # Calculate metabolic scores
  spm <- calculateMetabolicScores(
    spm,
    pathways = c("oxphos", "glycolysis"),
    verbose = FALSE
  )

  # Detect spatial patterns
  spm <- detectMetabolicGradients(spm, method = "moran", verbose = FALSE)

  # Compare conditions
  de_results <- compareMetabolicStates(
    spm,
    test_type = "t.test",
    verbose = FALSE
  )

  # Check that all components are present
  expect_true(nrow(metabolicScores(spm)) > 0)
  expect_true(length(analysisResults(spm)) > 0)
  expect_true(nrow(de_results) > 0)

  # Ensure object is still valid
  expect_true(validObject(spm))
})
