#' Detect spatial metabolic gradients
#'
#' @param object SpatialMetabolic object
#' @param features Features to test (default: all metabolic pathways)
#' @param method Method for spatial detection ("moran", "geary", "SPARK", "variogram")
#' @param n_neighbors Number of neighbors for spatial graph
#' @param graph_type Type of spatial graph ("knn", "delaunay", "radius")
#' @param radius Radius for radius-based graph (only used if graph_type = "radius")
#' @param permutations Number of permutations for significance testing
#' @param BPPARAM BiocParallel parameters
#'
#' @return SpatialMetabolic object with gradient results
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom stats dist pnorm p.adjust
#' @export
#'
#' @examples
#' \dontrun{
#' # Detect gradients in metabolic scores
#' spm <- detectMetabolicGradients(spm)
#'
#' # Test specific features
#' spm <- detectMetabolicGradients(spm,
#'                                features = c("Ndufa1", "Atp5a1"),
#'                                method = "moran")
#'
#' # Use different spatial graph
#' spm <- detectMetabolicGradients(spm,
#'                                method = "moran",
#'                                graph_type = "delaunay")
#' }
setMethod("detectMetabolicGradients", "SpatialMetabolic",
          function(object,
                   features = NULL,
                   method = "moran",
                   n_neighbors = 6,
                   graph_type = "knn",
                   radius = NULL,
                   permutations = 99,
                   BPPARAM = SerialParam()) {

            # Determine features to test
            if (is.null(features)) {
              if (nrow(metabolicScores(object)) > 0) {
                test_mat <- metabolicScores(object)
                feature_type <- "scores"
                message("Testing metabolic pathway scores for spatial patterns")
              } else {
                # Use highly variable genes if available
                if ("variable_feature" %in% colnames(rowData(object))) {
                  features <- rownames(object)[rowData(object)$variable_feature]
                  test_mat <- logcounts(object)[features, , drop = FALSE]
                  feature_type <- "genes"
                  message("Testing ", length(features), " variable genes")
                } else {
                  stop("No features specified and no metabolic scores calculated")
                }
              }
            } else {
              # Test specified features
              if (is.character(features)) {
                # Check if features are genes or pathways
                if (all(features %in% rownames(object))) {
                  test_mat <- logcounts(object)[features, , drop = FALSE]
                  feature_type <- "genes"
                } else if (all(features %in% rownames(metabolicScores(object)))) {
                  test_mat <- metabolicScores(object)[features, , drop = FALSE]
                  feature_type = "scores"
                } else {
                  stop("Features not found in genes or metabolic scores")
                }
              } else {
                stop("Features must be character vector")
              }
            }

            message("Testing ", nrow(test_mat), " ", feature_type, " for spatial patterns...")

            # Build spatial graph if not present
            graph_name <- paste0(graph_type, "_k", n_neighbors)
            if (!graph_name %in% names(spatialGraphs(object))) {
              message("Building ", graph_type, " spatial graph with k=", n_neighbors, "...")
              object <- .buildSpatialGraph(object, n_neighbors, graph_type, radius)
            }

            # Get the appropriate graph
            adj_mat <- spatialGraphs(object)[[graph_name]]

            # Calculate spatial statistics
            if (method == "moran") {
              results <- .calculateMoransI(test_mat, adj_mat, permutations, BPPARAM)
            } else if (method == "geary") {
              results <- .calculateGearysC(test_mat, adj_mat, permutations, BPPARAM)
            } else if (method == "variogram") {
              results <- .calculateVariogram(object, test_mat, BPPARAM)
            } else if (method == "SPARK") {
              results <- .runSPARK(object, test_mat, BPPARAM)
            } else {
              stop("Method must be one of: moran, geary, variogram, SPARK")
            }

            # Store results
            if (is.null(analysisResults(object))) {
              analysisResults(object) <- list()
            }

            result_name <- paste0("spatial_gradients_", feature_type)
            current_results <- analysisResults(object)
            current_results[[result_name]] <- results
            analysisResults(object) <- current_results

            # Identify significant features
            sig_features <- results$feature[results$adj_pvalue < 0.05]
            message("Found ", length(sig_features), " features with significant spatial patterns (FDR < 0.05)")

            if (length(sig_features) > 0 && length(sig_features) <= 10) {
              message("  Top features: ", paste(head(sig_features, 5), collapse = ", "))
            }

            return(object)
          })

#' Build spatial neighborhood graph
#' @keywords internal
#' @importFrom stats dist
.buildSpatialGraph <- function(object, n_neighbors, graph_type, radius) {
  coords <- spatialCoords(object)
  n_spots <- nrow(coords)

  if (graph_type == "knn") {
    # K-nearest neighbors graph
    dist_mat <- as.matrix(dist(coords))

    # Set diagonal to Inf to exclude self
    diag(dist_mat) <- Inf

    # Find k nearest neighbors for each spot
    nn_mat <- t(apply(dist_mat, 1, function(x) {
      order(x)[1:min(n_neighbors, n_spots - 1)]
    }))

    # Create adjacency matrix
    adj_mat <- matrix(0, nrow = n_spots, ncol = n_spots)
    for (i in seq_len(nrow(nn_mat))) {
      adj_mat[i, nn_mat[i, ]] <- 1
    }

    # Make symmetric (mutual nearest neighbors)
    adj_mat <- pmax(adj_mat, t(adj_mat))

  } else if (graph_type == "delaunay") {
    # Delaunay triangulation
    # Simplified implementation - for full implementation would use deldir or similar
    stop("Delaunay triangulation not yet implemented. Please use 'knn' or 'radius'.")

  } else if (graph_type == "radius") {
    # Radius-based graph
    if (is.null(radius)) {
      # Estimate radius as median k-NN distance
      dist_mat <- as.matrix(dist(coords))
      diag(dist_mat) <- Inf
      knn_dists <- apply(dist_mat, 1, function(x) sort(x)[n_neighbors])
      radius <- median(knn_dists)
      message("  Using radius = ", round(radius, 2))
    }

    dist_mat <- as.matrix(dist(coords))
    adj_mat <- (dist_mat <= radius & dist_mat > 0) * 1

  } else {
    stop("graph_type must be one of: knn, delaunay, radius")
  }

  # Store graph
  graph_name <- paste0(graph_type, "_k", n_neighbors)
  current_graphs <- spatialGraphs(object)
  current_graphs[[graph_name]] <- adj_mat
  spatialGraphs(object) <- current_graphs

  # Summary statistics
  n_edges <- sum(adj_mat) / 2  # Divided by 2 because symmetric
  avg_degree <- mean(rowSums(adj_mat))
  message("  Created graph with ", n_edges, " edges (avg degree: ",
          round(avg_degree, 1), ")")

  return(object)
}

#' Calculate Moran's I with permutation testing
#' @keywords internal
#' @importFrom stats var
.calculateMoransI <- function(test_mat, adj_mat, permutations, BPPARAM) {

  n <- ncol(test_mat)
  W <- sum(adj_mat)

  # Row normalize the weight matrix
  row_sums <- rowSums(adj_mat)
  row_sums[row_sums == 0] <- 1  # Avoid division by zero
  W_norm <- adj_mat / row_sums

  results <- bplapply(seq_len(nrow(test_mat)), function(i) {
    x <- as.numeric(test_mat[i, ])

    # Skip if no variance
    if (var(x) == 0) {
      return(data.frame(
        feature = rownames(test_mat)[i],
        morans_i = NA,
        expected_i = -1/(n-1),
        variance = NA,
        z_score = NA,
        pvalue = 1,
        pvalue_perm = 1,
        stringsAsFactors = FALSE
      ))
    }

    # Calculate observed Moran's I
    x_centered <- x - mean(x)

    # Numerator: sum of spatial covariance
    num <- sum(W_norm * outer(x_centered, x_centered))

    # Denominator: variance
    denom <- sum(x_centered^2) / n

    I_obs <- num / denom

    # Expected value under null hypothesis
    EI <- -1 / (n - 1)

    # Analytical variance (simplified)
    b2 <- n * sum(x_centered^4) / sum(x_centered^2)^2

    S1 <- 2 * sum(W_norm^2)
    S2 <- sum(rowSums(W_norm)^2)

    VI <- (n * ((n^2 - 3*n + 3) * S1 - n * S2 + 3 * W^2) -
             b2 * ((n^2 - n) * S1 - 2 * n * S2 + 6 * W^2)) /
      ((n - 1) * (n - 2) * (n - 3) * W^2)

    # Z-score and p-value
    z <- (I_obs - EI) / sqrt(VI)
    p_analytical <- 2 * pnorm(-abs(z))

    # Permutation test if requested
    if (permutations > 0) {
      I_perm <- replicate(permutations, {
        x_perm <- sample(x)
        x_perm_centered <- x_perm - mean(x_perm)
        sum(W_norm * outer(x_perm_centered, x_perm_centered)) /
          (sum(x_perm_centered^2) / n)
      })

      # Two-tailed p-value
      p_perm <- (sum(abs(I_perm - EI) >= abs(I_obs - EI)) + 1) / (permutations + 1)
    } else {
      p_perm <- p_analytical
    }

    return(data.frame(
      feature = rownames(test_mat)[i],
      morans_i = I_obs,
      expected_i = EI,
      variance = VI,
      z_score = z,
      pvalue = p_analytical,
      pvalue_perm = p_perm,
      stringsAsFactors = FALSE
    ))
  }, BPPARAM = BPPARAM)

  results <- do.call(rbind, results)

  # Use permutation p-values if available
  if (permutations > 0) {
    results$pvalue <- results$pvalue_perm
  }

  results$adj_pvalue <- p.adjust(results$pvalue, method = "BH")

  # Sort by significance
  results <- results[order(results$adj_pvalue, -abs(results$morans_i)), ]

  return(results)
}

#' Calculate Geary's C
#' @keywords internal
.calculateGearysC <- function(test_mat, adj_mat, permutations, BPPARAM) {

  n <- ncol(test_mat)
  W <- sum(adj_mat)

  results <- bplapply(seq_len(nrow(test_mat)), function(i) {
    x <- as.numeric(test_mat[i, ])

    # Skip if no variance
    if (var(x) == 0) {
      return(data.frame(
        feature = rownames(test_mat)[i],
        gearys_c = NA,
        expected_c = 1,
        z_score = NA,
        pvalue = 1,
        stringsAsFactors = FALSE
      ))
    }

    # Calculate Geary's C
    x_mean <- mean(x)

    # Numerator: sum of squared differences between neighbors
    num <- 0
    for (j in 1:n) {
      for (k in 1:n) {
        if (adj_mat[j, k] == 1) {
          num <- num + (x[j] - x[k])^2
        }
      }
    }
    num <- num * (n - 1)

    # Denominator: total sum of squares
    denom <- 2 * W * sum((x - x_mean)^2)

    C <- num / denom

    # Expected value
    EC <- 1

    # Simplified z-score (approximation)
    # For exact variance calculation, see Cliff & Ord (1981)
    z <- (C - EC) / sqrt(2 / (n * W))
    p <- 2 * pnorm(-abs(z))

    return(data.frame(
      feature = rownames(test_mat)[i],
      gearys_c = C,
      expected_c = EC,
      z_score = z,
      pvalue = p,
      stringsAsFactors = FALSE
    ))
  }, BPPARAM = BPPARAM)

  results <- do.call(rbind, results)
  results$adj_pvalue <- p.adjust(results$pvalue, method = "BH")
  results <- results[order(results$adj_pvalue), ]

  return(results)
}

#' Calculate variogram-based spatial statistics
#' @keywords internal
.calculateVariogram <- function(object, test_mat, BPPARAM) {
  coords <- spatialCoords(object)

  # Calculate pairwise distances
  dist_mat <- as.matrix(dist(coords))

  # Define distance bins
  max_dist <- max(dist_mat[upper.tri(dist_mat)])
  n_bins <- 20
  breaks <- seq(0, max_dist, length.out = n_bins + 1)

  results <- bplapply(seq_len(nrow(test_mat)), function(i) {
    x <- as.numeric(test_mat[i, ])

    if (var(x) == 0) {
      return(data.frame(
        feature = rownames(test_mat)[i],
        range = NA,
        sill = NA,
        nugget = NA,
        spatial_ratio = NA,
        stringsAsFactors = FALSE
      ))
    }

    # Calculate empirical variogram
    n <- length(x)
    gamma <- numeric(n_bins)
    counts <- numeric(n_bins)

    for (j in 1:(n-1)) {
      for (k in (j+1):n) {
        dist_jk <- dist_mat[j, k]
        bin <- which(breaks[-1] >= dist_jk)[1]
        if (!is.na(bin) && bin <= n_bins) {
          gamma[bin] <- gamma[bin] + 0.5 * (x[j] - x[k])^2
          counts[bin] <- counts[bin] + 1
        }
      }
    }

    # Average by number of pairs
    gamma[counts > 0] <- gamma[counts > 0] / counts[counts > 0]

    # Simple variogram parameters
    # Range: distance at which variance plateaus (simplified)
    mid_bins <- (breaks[-1] + breaks[-length(breaks)]) / 2
    valid_bins <- which(counts > 10)  # Only use bins with enough pairs

    if (length(valid_bins) >= 3) {
      # Estimate sill as median of last few bins
      sill <- median(gamma[tail(valid_bins, 3)])

      # Find range as distance where gamma reaches 95% of sill
      range_idx <- which(gamma[valid_bins] >= 0.95 * sill)[1]
      if (!is.na(range_idx)) {
        range <- mid_bins[valid_bins[range_idx]]
      } else {
        range <- max_dist / 2
      }

      # Nugget as y-intercept (first bin value)
      nugget <- gamma[valid_bins[1]]

      # Spatial ratio: proportion of variance that is spatially structured
      spatial_ratio <- (sill - nugget) / sill
    } else {
      range <- NA
      sill <- var(x)
      nugget <- 0
      spatial_ratio <- NA
    }

    return(data.frame(
      feature = rownames(test_mat)[i],
      range = range,
      sill = sill,
      nugget = nugget,
      spatial_ratio = spatial_ratio,
      stringsAsFactors = FALSE
    ))
  }, BPPARAM = BPPARAM)

  results <- do.call(rbind, results)

  # Sort by spatial ratio (higher = more spatial structure)
  results <- results[order(-results$spatial_ratio, na.last = TRUE), ]

  return(results)
}

#' Placeholder for SPARK integration
#' @keywords internal
.runSPARK <- function(object, test_mat, BPPARAM) {
  stop("SPARK integration not yet implemented. Please use 'moran', 'geary', or 'variogram'.")
}

#' Find metabolic neighborhoods with distinct signatures
#'
#' @param object SpatialMetabolic object
#' @param k Number of neighborhoods to identify
#' @param method Clustering method ("kmeans", "leiden", "louvain")
#' @param use_pca Whether to use PCA before clustering
#' @param n_pcs Number of PCs to use
#' @param resolution Resolution parameter for graph-based clustering
#' @param seed Random seed
#'
#' @return SpatialMetabolic object with neighborhood assignments
#' @importFrom stats kmeans prcomp
#' @export
setMethod("findMetabolicNeighborhoods", "SpatialMetabolic",
          function(object,
                   k = 5,
                   method = "kmeans",
                   use_pca = TRUE,
                   n_pcs = 10,
                   resolution = 0.8,
                   seed = 42) {

            # Check for metabolic scores
            if (nrow(metabolicScores(object)) == 0) {
              stop("No metabolic scores found. Run calculateMetabolicScores first.")
            }

            # Get data for clustering
            cluster_data <- t(metabolicScores(object))

            # PCA if requested
            if (use_pca) {
              if (n_pcs > ncol(cluster_data)) {
                n_pcs <- ncol(cluster_data)
              }

              pca <- prcomp(cluster_data, scale. = TRUE, center = TRUE)
              cluster_data <- pca$x[, 1:n_pcs]
            }

            set.seed(seed)

            # Perform clustering
            if (method == "kmeans") {
              km <- kmeans(cluster_data, centers = k, nstart = 50)
              neighborhoods <- factor(km$cluster)
            } else if (method %in% c("leiden", "louvain")) {
              stop("Graph-based clustering not yet implemented. Please use 'kmeans'.")
            } else {
              stop("Method must be one of: kmeans, leiden, louvain")
            }

            # Add to colData
            colData(object)$metabolic_neighborhood <- neighborhoods

            # Store in analysis results
            current_results <- analysisResults(object)
            current_results[["metabolic_neighborhoods"]] <- list(
              assignments = neighborhoods,
              method = method,
              k = k,
              centers = if (method == "kmeans") km$centers else NULL
            )
            analysisResults(object) <- current_results

            message("Identified ", k, " metabolic neighborhoods using ", method)

            # Summary
            table_neighborhoods <- table(neighborhoods)
            for (i in seq_len(k)) {
              message("  Neighborhood ", i, ": ", table_neighborhoods[i], " spots")
            }

            return(object)
          })
