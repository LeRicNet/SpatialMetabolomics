# tests/testthat/helper-test-data.R

# Helper function to create test data
create_test_data <- function(n_genes = 100, n_spots = 50, n_samples = 2) {
  # Handle edge case of zero dimensions
  if (n_genes == 0 || n_spots == 0) {
    # Create truly empty matrices and data structures
    counts <- matrix(integer(0), nrow = n_genes, ncol = n_spots * n_samples)

    # Set appropriate dimnames
    if (n_genes > 0) {
      rownames(counts) <- paste0("Gene", seq_len(n_genes))
    }
    if (n_spots > 0) {
      colnames(counts) <- unlist(lapply(seq_len(n_samples), function(i) {
        paste0("Spot", seq_len(n_spots), "_S", i)
      }))
    }

    # Create empty metadata
    col_data <- data.frame(
      sample_id = character(0),
      condition = character(0),
      row.names = character(0)
    )

    # Create empty coordinates
    coords <- matrix(numeric(0), nrow = n_spots * n_samples, ncol = 2)
    colnames(coords) <- c("x", "y")

  } else {
    # Normal case with non-zero dimensions
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

    # Create spatial coordinates - fix to respect n_spots parameter
    coords <- do.call(rbind, lapply(seq_len(n_samples), function(i) {
      # Create a grid pattern that respects n_spots
      if (n_spots <= 50) {
        # For small n_spots, create a simple grid
        grid_size <- ceiling(sqrt(n_spots))
        x_coords <- rep(1:grid_size, length.out = n_spots)
        y_coords <- rep(1:grid_size, each = grid_size)[1:n_spots]
      } else {
        # For larger n_spots, use the original 10x5 pattern repeated
        base_x <- rep(1:10, 5)
        base_y <- rep(1:5, each = 10)
        x_coords <- rep(base_x, length.out = n_spots)
        y_coords <- rep(base_y, length.out = n_spots)
      }

      data.frame(
        x = x_coords,
        y = y_coords
      )
    }))
  }

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

# Helper to check if plot is valid ggplot
expect_ggplot <- function(p) {
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_true(length(p$layers) > 0)
  testthat::expect_true(!is.null(p$data) || !is.null(p$layers[[1]]$data))
}
