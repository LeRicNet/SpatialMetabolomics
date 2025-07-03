# tests/testthat/helper-test-data.R

# Helper function to create test data
create_test_data <- function(n_genes = 100, n_spots = 50, n_samples = 2) {
  # Handle edge cases
  if (n_genes == 0 || n_spots == 0) {
    # Create minimal valid object
    if (n_genes == 0) n_genes <- 1
    if (n_spots == 0) n_spots <- 1
  }

  # Create gene names that will work with test pathways
  all_gene_names <- paste0("Gene", seq_len(n_genes))

  counts_list <- lapply(seq_len(n_samples), function(i) {
    matrix(rpois(n_genes * n_spots, lambda = 5), nrow = n_genes, ncol = n_spots,
           dimnames = list(all_gene_names, paste0("Spot", seq_len(n_spots), "_S", i)))
  })

  counts <- do.call(cbind, counts_list)
  col_data <- data.frame(
    sample_id = rep(paste0("Sample", seq_len(n_samples)), each = n_spots),
    condition = rep(c("WT", "Disease"), length.out = n_spots * n_samples),
    row.names = colnames(counts)
  )

  coords <- do.call(rbind, lapply(seq_len(n_samples), function(i) {
    data.frame(x = rep(1:10, length.out = n_spots), y = rep(1:5, each = 10)[seq_len(n_spots)])
  }))

  spe <- SpatialExperiment(assays = list(counts = counts), colData = col_data, spatialCoords = as.matrix(coords))
  spm <- as(spe, "SpatialMetabolic")
  comparisonGroups(spm) <- factor(col_data$condition, levels = c("WT", "Disease"))
  return(spm)
}

# Helper to check if plot is valid ggplot
expect_ggplot <- function(p) {
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_true(length(p$layers) > 0)
  testthat::expect_true(!is.null(p$data) || !is.null(p$layers[[1]]$data))
}
