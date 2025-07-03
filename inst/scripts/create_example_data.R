#' Create example data for SpatialMetabolics package
#'
#' This script creates a small example dataset for package examples and testing
#'
#' @author Your Name

library(SpatialExperiment)
library(SingleCellExperiment)

# Set seed for reproducibility
set.seed(42)

# Parameters
n_genes <- 1000
n_spots_per_sample <- 500  # Small for examples
n_samples <- 4

# Create gene names including real metabolic genes
gene_names <- c(
  # Include some real OXPHOS genes
  "Ndufa1", "Ndufa2", "Ndufa3", "Ndufa4", "Ndufa5",
  "Sdha", "Sdhb", "Sdhc", "Sdhd",
  "Uqcrc1", "Uqcrc2", "Uqcrfs1",
  "Cox4i1", "Cox5a", "Cox6a1", "Cox7a1",
  "Atp5a1", "Atp5b", "Atp5c1", "Atp5d",

  # Glycolysis genes
  "Hk1", "Hk2", "Pfkl", "Pfkm", "Aldoa",
  "Gapdh", "Pgk1", "Eno1", "Pkm", "Ldha",

  # Mitochondrial biogenesis
  "Ppargc1a", "Nrf1", "Tfam", "Polrmt", "Sirt3",

  # Random genes
  paste0("Gene", 36:n_genes)
)

# Function to create sample data
create_sample <- function(sample_id, condition) {
  # Base counts - Poisson distributed
  base_lambda <- runif(n_genes, min = 1, max = 10)
  counts <- matrix(
    rpois(n_genes * n_spots_per_sample, lambda = rep(base_lambda, n_spots_per_sample)),
    nrow = n_genes,
    ncol = n_spots_per_sample
  )

  # Add condition-specific effects
  if (condition == "Disease") {
    # Decrease OXPHOS genes (indices 1-20)
    counts[1:20, ] <- round(counts[1:20, ] * runif(20, 0.5, 0.8))

    # Increase glycolysis genes (indices 21-30)
    counts[21:30, ] <- round(counts[21:30, ] * runif(10, 1.2, 1.8))

    # Add some spatial pattern
    x_coords <- rep(1:25, 20)[1:n_spots_per_sample]
    y_coords <- rep(1:20, each = 25)[1:n_spots_per_sample]

    # Create gradient for some genes
    for (i in 1:5) {
      spatial_effect <- 1 + 0.5 * sin(x_coords/5) * cos(y_coords/5)
      counts[i, ] <- round(counts[i, ] * spatial_effect)
    }
  }

  # Add technical noise
  counts <- counts + matrix(
    rpois(n_genes * n_spots_per_sample, lambda = 0.1),
    nrow = n_genes,
    ncol = n_spots_per_sample
  )

  # Set gene and spot names
  rownames(counts) <- gene_names[1:n_genes]
  colnames(counts) <- paste0(sample_id, "_Spot", 1:n_spots_per_sample)

  # Create spatial coordinates (grid layout)
  coords <- data.frame(
    x = rep(1:25, 20)[1:n_spots_per_sample],
    y = rep(1:20, each = 25)[1:n_spots_per_sample]
  )

  # Create metadata
  metadata <- data.frame(
    sample_id = sample_id,
    condition = condition,
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  return(list(
    counts = counts,
    coords = coords,
    metadata = metadata
  ))
}

# Create data for each sample
samples <- list(
  WT1 = create_sample("WT1", "Control"),
  WT2 = create_sample("WT2", "Control"),
  Disease1 = create_sample("Disease1", "Disease"),
  Disease2 = create_sample("Disease2", "Disease")
)

# Combine all data
all_counts <- do.call(cbind, lapply(samples, "[[", "counts"))
all_coords <- do.call(rbind, lapply(samples, "[[", "coords"))
all_metadata <- do.call(rbind, lapply(samples, "[[", "metadata"))

# Create SpatialExperiment object
spe <- SpatialExperiment(
  assays = list(counts = as(all_counts, "dgCMatrix")),  # Sparse matrix
  colData = all_metadata,
  spatialCoords = as.matrix(all_coords)
)

# Add some additional metadata
colData(spe)$total_counts <- colSums(counts(spe))
colData(spe)$n_genes <- colSums(counts(spe) > 0)

# Calculate simple QC metrics
mt_genes <- grep("^mt-|^Mt-", rownames(spe), value = TRUE)
if (length(mt_genes) > 0) {
  colData(spe)$percent_mt <- colSums(counts(spe)[mt_genes, ]) / colData(spe)$total_counts * 100
} else {
  colData(spe)$percent_mt <- 0
}

# Save as RDS file
saveRDS(spe, "spatialmetabolics_example_data.rds")

# Create a smaller subset for quick examples
spe_small <- spe[1:500, seq(1, ncol(spe), by = 2)]
saveRDS(spe_small, "spatialmetabolics_example_data_small.rds")

# Print summary
cat("Created example SpatialExperiment object with:\n")
cat(" -", nrow(spe), "genes\n")
cat(" -", ncol(spe), "spots\n")
cat(" -", length(unique(colData(spe)$sample_id)), "samples\n")
cat(" -", length(unique(colData(spe)$condition)), "conditions\n")
cat("\nFiles saved:\n")
cat(" - spatialmetabolics_example_data.rds (full)\n")
cat(" - spatialmetabolics_example_data_small.rds (subset)\n")
