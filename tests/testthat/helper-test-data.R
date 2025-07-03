# tests/testthat/helper-test-data.R

# Helper function to create test data for all tests
create_test_data <- function(n_genes = 100, n_spots = 50, n_samples = 2) {

  # Handle edge cases
  if (n_genes == 0 || n_spots == 0 || n_samples == 0) {
    # Return minimal valid object
    counts <- matrix(integer(0), nrow = 0, ncol = 0)
    col_data <- data.frame(
      sample_id = character(0),
      condition = character(0)
    )
    coords <- matrix(numeric(0), nrow = 0, ncol = 2)
    colnames(coords) <- c("x", "y")

    spe <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = counts),
      colData = col_data,
      spatialCoords = coords
    )

    spm <- methods::as(spe, "SpatialMetabolic")
    comparisonGroups(spm) <- factor(col_data$condition, levels = c("WT", "Disease"))
    return(spm)
  }

  # SPECIAL CASE: For single gene/single spot tests, use only 1 sample
  if (n_genes == 1 && n_spots == 1) {
    n_samples <- 1
  }

  # Create count matrices with some real pathway genes
  real_genes <- c(
    "Ndufa1", "Ndufa2", "Atp5a1", "Atp5b", "Cox4i1", "Cox5a",  # OXPHOS
    "Hk1", "Hk2", "Pfkl", "Gapdh", "Pgk1", "Pkm", "Ldha",      # Glycolysis
    paste0("TestGene", 1:max(0, n_genes - 13))                   # Fill remainder
  )

  gene_names <- real_genes[1:n_genes]

  counts_list <- lapply(seq_len(n_samples), function(i) {
    matrix(
      rpois(n_genes * n_spots, lambda = 5),
      nrow = n_genes,
      ncol = n_spots,
      dimnames = list(
        gene_names,
        paste0("Spot", seq_len(n_spots), "_S", i)
      )
    )
  })

  # Combine counts
  counts <- do.call(cbind, counts_list)

  # Create metadata with explicit factor levels
  total_spots <- n_spots * n_samples

  # Handle condition assignment for different sample numbers
  if (n_samples == 1) {
    conditions <- rep("WT", total_spots)
  } else {
    conditions <- rep(c("WT", "Disease"), length.out = total_spots)
  }

  col_data <- data.frame(
    sample_id = rep(paste0("Sample", seq_len(n_samples)), each = n_spots),
    condition = conditions,
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  # Ensure factor has correct levels
  col_data$condition <- factor(col_data$condition, levels = c("WT", "Disease"))

  # Create spatial coordinates - fix dimension matching
  coords_list <- lapply(seq_len(n_samples), function(i) {
    # Create a grid that can accommodate n_spots
    if (n_spots == 1) {
      grid <- data.frame(x = 1, y = 1)
    } else {
      grid_size <- ceiling(sqrt(n_spots))
      grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)[1:n_spots, ]
    }
    return(grid)
  })

  coords <- do.call(rbind, coords_list)
  rownames(coords) <- colnames(counts)

  # Create SpatialExperiment
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = col_data,
    spatialCoords = as.matrix(coords)
  )

  # Convert to SpatialMetabolic
  spm <- methods::as(spe, "SpatialMetabolic")
  comparisonGroups(spm) <- col_data$condition  # Already has correct levels

  return(spm)
}

# Helper to check if plot is valid ggplot
expect_ggplot <- function(p) {
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_true(length(p$layers) > 0)
  testthat::expect_true(!is.null(p$data) || !is.null(p$layers[[1]]$data))
}
