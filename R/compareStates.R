#' Compare metabolic states between conditions
#'
#' @param object SpatialMetabolic object
#' @param condition Column name with condition information
#' @param ref_group Reference group name
#' @param test_group Test group name
#' @param test_type Type of test ("wilcox", "t.test", "spatial_lm", "permutation")
#' @param features Features to test (default: metabolic scores)
#' @param min_spots Minimum spots per group for testing
#' @param spatial_correction Whether to correct for spatial autocorrelation
#' @param n_permutations Number of permutations for permutation test
#' @param verbose Print progress messages
#' @param BPPARAM BiocParallel parameters
#'
#' @return Data frame with differential expression results
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom stats wilcox.test t.test p.adjust sd var
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment colData
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare metabolic scores between conditions
#' results <- compareMetabolicStates(spm,
#'                                  condition = "condition",
#'                                  ref_group = "WT",
#'                                  test_group = "SVD")
#'
#' # Compare specific genes
#' gene_results <- compareMetabolicStates(spm,
#'                                       condition = "condition",
#'                                       ref_group = "WT",
#'                                       test_group = "SVD",
#'                                       features = c("Ndufa1", "Atp5a1"),
#'                                       test_type = "t.test")
#' }
setMethod("compareMetabolicStates", "SpatialMetabolic",
          function(object,
                   condition = "condition",
                   ref_group = NULL,
                   test_group = NULL,
                   test_type = "wilcox",
                   features = NULL,
                   min_spots = 3,
                   spatial_correction = FALSE,
                   n_permutations = 1000,
                   verbose = TRUE,
                   BPPARAM = SerialParam()) {

            # Get condition information
            conditions <- colData(object)[[condition]]

            if (is.null(conditions)) {
              stop("Condition column '", condition, "' not found in colData")
            }

            # Convert to factor if needed
            if (!is.factor(conditions)) {
              conditions <- factor(conditions)
            }

            # Determine groups
            unique_conditions <- levels(conditions)

            if (is.null(ref_group) || is.null(test_group)) {
              if (length(unique_conditions) != 2) {
                stop("Must specify ref_group and test_group when more than 2 conditions exist")
              }
              ref_group <- unique_conditions[1]
              test_group <- unique_conditions[2]
              if (verbose) {
                message("Comparing ", test_group, " vs ", ref_group)
              }
            }

            # Validate groups exist
            if (!ref_group %in% unique_conditions) {
              stop("Reference group '", ref_group, "' not found in conditions")
            }
            if (!test_group %in% unique_conditions) {
              stop("Test group '", test_group, "' not found in conditions")
            }

            # Subset to groups of interest
            keep_cells <- conditions %in% c(ref_group, test_group)
            object_sub <- object[, keep_cells]
            conditions_sub <- droplevels(conditions[keep_cells])

            # Check group sizes
            group_sizes <- table(conditions_sub)
            if (any(group_sizes < min_spots)) {
              stop("Groups have too few spots: ",
                   paste(names(group_sizes), "=", group_sizes, collapse = ", "))
            }

            # Determine features to test
            if (is.null(features)) {
              if (nrow(metabolicScores(object)) > 0) {
                test_mat <- t(metabolicScores(object_sub))
                feature_type <- "pathways"
                feature_names <- colnames(test_mat)
              } else {
                # Test all genes
                test_mat <- t(logcounts(object_sub))
                feature_type <- "genes"
                feature_names <- colnames(test_mat)

                # Limit to expressed genes
                expression_rate <- colMeans(test_mat > 0)
                expressed_genes <- names(expression_rate)[expression_rate > 0.1]
                test_mat <- test_mat[, expressed_genes, drop = FALSE]
                feature_names <- colnames(test_mat)

                if (verbose) {
                  message("Testing ", length(feature_names), " expressed genes")
                }
              }
            } else {
              # Test specified features
              if (all(features %in% rownames(metabolicScores(object_sub)))) {
                score_mat <- metabolicScores(object_sub)[features, , drop = FALSE]
                # Ensure it's a matrix even with one feature
                if (!is.matrix(score_mat)) {
                  score_mat <- matrix(score_mat, nrow = 1,
                                      dimnames = list(features, colnames(metabolicScores(object_sub))))
                }
                test_mat <- t(score_mat)
                feature_type <- "pathways"
              } else if (all(features %in% rownames(object_sub))) {
                expr_mat <- logcounts(object_sub)[features, , drop = FALSE]
                # Ensure it's a matrix even with one feature
                if (!is.matrix(expr_mat)) {
                  expr_mat <- matrix(expr_mat, nrow = 1,
                                     dimnames = list(features, colnames(logcounts(object_sub))))
                }
                test_mat <- t(expr_mat)
                feature_type <- "genes"
              } else {
                # Check which features are missing from the subsetted object
                missing_scores <- setdiff(features, rownames(metabolicScores(object_sub)))
                missing_genes <- setdiff(features, rownames(object_sub))

                if (length(missing_scores) < length(missing_genes)) {
                  # More features available in metabolic scores
                  available_features <- intersect(features, rownames(metabolicScores(object_sub)))
                  if (length(available_features) == 0) {
                    stop("No features found in metabolic scores")
                  }
                  if (verbose) {
                    message("Using ", length(available_features), " available features from metabolic scores")
                  }
                  score_mat <- metabolicScores(object_sub)[available_features, , drop = FALSE]
                  if (!is.matrix(score_mat)) {
                    score_mat <- matrix(score_mat, nrow = 1,
                                        dimnames = list(available_features, colnames(metabolicScores(object_sub))))
                  }
                  test_mat <- t(score_mat)
                  feature_type <- "pathways"
                  features <- available_features
                } else {
                  # More features available in gene expression
                  available_features <- intersect(features, rownames(object_sub))
                  if (length(available_features) == 0) {
                    stop("No features found in gene expression data")
                  }
                  if (verbose) {
                    message("Using ", length(available_features), " available features from gene expression")
                  }
                  expr_mat <- logcounts(object_sub)[available_features, , drop = FALSE]
                  if (!is.matrix(expr_mat)) {
                    expr_mat <- matrix(expr_mat, nrow = 1,
                                       dimnames = list(available_features, colnames(logcounts(object_sub))))
                  }
                  test_mat <- t(expr_mat)
                  feature_type <- "genes"
                  features <- available_features
                }
              }
              feature_names <- colnames(test_mat)
            }

            if (verbose) {
              message("Testing ", ncol(test_mat), " ", feature_type, " using ", test_type, " test...")
            }

            # Get indices for each group
            ref_idx <- which(conditions_sub == ref_group)
            test_idx <- which(conditions_sub == test_group)

            # Perform tests
            if (test_type == "wilcox") {
              results <- .wilcoxonTest(test_mat, ref_idx, test_idx, feature_names, BPPARAM)
            } else if (test_type == "t.test") {
              results <- .tTest(test_mat, ref_idx, test_idx, feature_names, BPPARAM)
            } else if (test_type == "permutation") {
              results <- .permutationTest(test_mat, ref_idx, test_idx, feature_names,
                                          n_permutations, BPPARAM)
            } else if (test_type == "spatial_lm") {
              if (!spatial_correction) {
                warning("spatial_lm requires spatial_correction=TRUE. Setting it now.")
                spatial_correction <- TRUE
              }
              results <- .spatialLinearModel(object_sub, test_mat, conditions_sub,
                                             ref_group, test_group, feature_names, BPPARAM)
            } else {
              stop("test_type must be one of: wilcox, t.test, permutation, spatial_lm")
            }

            # Add fold change and mean expression
            ref_means <- colMeans(test_mat[ref_idx, , drop = FALSE])
            test_means <- colMeans(test_mat[test_idx, , drop = FALSE])

            results$mean_ref <- ref_means[results$feature]
            results$mean_test <- test_means[results$feature]
            results$log2FC <- log2(test_means[results$feature] + 1) - log2(ref_means[results$feature] + 1)
            results$mean_expression <- (ref_means[results$feature] + test_means[results$feature]) / 2

            # Calculate effect size (Cohen's d)
            cohens_d <- sapply(results$feature, function(f) {
              ref_vals <- test_mat[ref_idx, f]
              test_vals <- test_mat[test_idx, f]

              pooled_sd <- sqrt(((length(ref_vals) - 1) * var(ref_vals) +
                                   (length(test_vals) - 1) * var(test_vals)) /
                                  (length(ref_vals) + length(test_vals) - 2))

              if (pooled_sd == 0) return(0)

              (mean(test_vals) - mean(ref_vals)) / pooled_sd
            })

            results$cohens_d <- cohens_d

            # Add group information
            results$ref_group <- ref_group
            results$test_group <- test_group
            results$comparison <- paste0(test_group, "_vs_", ref_group)

            # Reorder columns
            col_order <- c("feature", "comparison", "ref_group", "test_group",
                           "mean_ref", "mean_test", "log2FC", "mean_expression",
                           "statistic", "cohens_d", "pvalue", "adj_pvalue")

            results <- results[, intersect(col_order, colnames(results))]

            # Sort by p-value
            results <- results[order(results$pvalue), ]
            rownames(results) <- NULL

            # Store results
            current_results <- analysisResults(object)
            comparison_name <- paste0(test_group, "_vs_", ref_group, "_", feature_type)
            current_results[[comparison_name]] <- results
            analysisResults(object) <- current_results

            # Summary
            sig_up <- sum(results$adj_pvalue < 0.05 & results$log2FC > 0, na.rm = TRUE)
            sig_down <- sum(results$adj_pvalue < 0.05 & results$log2FC < 0, na.rm = TRUE)

            if (verbose) {
              message("Found ", sig_up, " upregulated and ", sig_down,
                      " downregulated ", feature_type, " (FDR < 0.05)")

              if (sig_up + sig_down > 0 && sig_up + sig_down <= 10) {
                top_features <- head(results$feature[results$adj_pvalue < 0.05], 5)
                message("  Top significant: ", paste(top_features, collapse = ", "))
              }
            }

            return(results)
          })

#' Wilcoxon rank-sum test implementation
#' @keywords internal
.wilcoxonTest <- function(test_mat, ref_idx, test_idx, feature_names, BPPARAM) {

  results <- bplapply(seq_len(ncol(test_mat)), function(i) {
    ref_vals <- test_mat[ref_idx, i]
    test_vals <- test_mat[test_idx, i]

    # Skip if no variance in both groups
    if (sd(ref_vals) == 0 && sd(test_vals) == 0) {
      return(data.frame(
        feature = feature_names[i],
        statistic = NA,
        pvalue = 1,
        method = "wilcox",
        stringsAsFactors = FALSE
      ))
    }

    # Wilcoxon test
    tryCatch({
      wt <- wilcox.test(test_vals, ref_vals, exact = FALSE)

      data.frame(
        feature = feature_names[i],
        statistic = wt$statistic,
        pvalue = wt$p.value,
        method = "wilcox",
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        feature = feature_names[i],
        statistic = NA,
        pvalue = 1,
        method = "wilcox",
        stringsAsFactors = FALSE
      )
    })
  }, BPPARAM = BPPARAM)

  results <- do.call(rbind, results)
  results$adj_pvalue <- p.adjust(results$pvalue, method = "BH")

  return(results)
}

#' t-test implementation
#' @keywords internal
.tTest <- function(test_mat, ref_idx, test_idx, feature_names, BPPARAM) {

  results <- bplapply(seq_len(ncol(test_mat)), function(i) {
    ref_vals <- test_mat[ref_idx, i]
    test_vals <- test_mat[test_idx, i]

    # Skip if no variance in both groups
    if (sd(ref_vals) == 0 && sd(test_vals) == 0) {
      return(data.frame(
        feature = feature_names[i],
        statistic = NA,
        pvalue = 1,
        method = "t.test",
        stringsAsFactors = FALSE
      ))
    }

    # t-test
    tryCatch({
      tt <- t.test(test_vals, ref_vals, var.equal = FALSE)

      data.frame(
        feature = feature_names[i],
        statistic = tt$statistic,
        pvalue = tt$p.value,
        method = "t.test",
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        feature = feature_names[i],
        statistic = NA,
        pvalue = 1,
        method = "t.test",
        stringsAsFactors = FALSE
      )
    })
  }, BPPARAM = BPPARAM)

  results <- do.call(rbind, results)
  results$adj_pvalue <- p.adjust(results$pvalue, method = "BH")

  return(results)
}

#' Permutation test implementation
#' @keywords internal
.permutationTest <- function(test_mat, ref_idx, test_idx, feature_names,
                             n_permutations, BPPARAM) {

  all_idx <- c(ref_idx, test_idx)
  n_ref <- length(ref_idx)
  n_test <- length(test_idx)

  results <- bplapply(seq_len(ncol(test_mat)), function(i) {
    vals <- test_mat[all_idx, i]

    # Observed difference
    obs_diff <- mean(vals[(n_ref + 1):(n_ref + n_test)]) - mean(vals[1:n_ref])

    # Skip if no variance
    if (sd(vals) == 0) {
      return(data.frame(
        feature = feature_names[i],
        statistic = obs_diff,
        pvalue = 1,
        method = "permutation",
        stringsAsFactors = FALSE
      ))
    }

    # Permutation distribution
    perm_diffs <- replicate(n_permutations, {
      perm_idx <- sample(length(vals))
      mean(vals[perm_idx[(n_ref + 1):(n_ref + n_test)]]) -
        mean(vals[perm_idx[1:n_ref]])
    })

    # Two-tailed p-value
    pvalue <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_permutations + 1)

    data.frame(
      feature = feature_names[i],
      statistic = obs_diff,
      pvalue = pvalue,
      method = "permutation",
      stringsAsFactors = FALSE
    )
  }, BPPARAM = BPPARAM)

  results <- do.call(rbind, results)
  results$adj_pvalue <- p.adjust(results$pvalue, method = "BH")

  return(results)
}

#' Spatial linear model (placeholder)
#' @keywords internal
.spatialLinearModel <- function(object, test_mat, conditions, ref_group, test_group,
                                feature_names, BPPARAM) {
  stop("Spatial linear models not yet implemented. Please use 'wilcox', 't.test', or 'permutation'.")
}

#' Perform pathway set enrichment analysis
#'
#' @param de_results Results from compareMetabolicStates
#' @param pathways List of pathway gene sets
#' @param universe Character vector of background genes
#' @param pvalue_cutoff P-value cutoff for DE genes
#' @param log2fc_cutoff Log2 fold change cutoff
#' @param min_genes Minimum genes in pathway
#' @param method Enrichment method ("hypergeometric", "gsea")
#'
#' @return Data frame with enrichment results
#' @importFrom stats phyper p.adjust
#' @export
#'
#' @examples
#' \dontrun{
#' # Run enrichment on DE results
#' de_genes <- compareMetabolicStates(spm, features = rownames(spm))
#' pathways <- getAllMetabolicPathways("mouse")
#' enrichment <- pathwayEnrichment(de_genes, pathways)
#' }
pathwayEnrichment <- function(de_results,
                              pathways,
                              universe = NULL,
                              pvalue_cutoff = 0.05,
                              log2fc_cutoff = 0,
                              min_genes = 3,
                              method = "hypergeometric") {

  if (!"feature" %in% colnames(de_results)) {
    stop("de_results must have 'feature' column")
  }

  # Get gene universe
  if (is.null(universe)) {
    universe <- unique(de_results$feature)
  }

  # Get significant genes
  sig_up <- de_results$feature[
    de_results$adj_pvalue < pvalue_cutoff &
      de_results$log2FC > log2fc_cutoff
  ]

  sig_down <- de_results$feature[
    de_results$adj_pvalue < pvalue_cutoff &
      de_results$log2FC < -log2fc_cutoff
  ]

  # FIX: Check if we have any significant genes
  if (length(sig_up) == 0 && length(sig_down) == 0) {
    message("No significant genes found with the given cutoffs")
    return(data.frame(
      pathway = character(0),
      direction = character(0),
      n_pathway_genes = integer(0),
      n_sig_genes = integer(0),
      n_overlap = integer(0),
      expected_overlap = numeric(0),
      fold_enrichment = numeric(0),
      pvalue = numeric(0),
      genes = character(0),
      adj_pvalue = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # Run enrichment for up and down separately
  if (method == "hypergeometric") {
    results_up <- .hypergeometricTest(sig_up, pathways, universe,
                                      min_genes, direction = "up")
    results_down <- .hypergeometricTest(sig_down, pathways, universe,
                                        min_genes, direction = "down")

    results <- rbind(results_up, results_down)
  } else if (method == "gsea") {
    stop("GSEA method not yet implemented")
  } else {
    stop("Method must be 'hypergeometric' or 'gsea'")
  }

  # FIX: Check if results is empty before sorting
  if (is.null(results) || nrow(results) == 0) {
    message("No enriched pathways found")
    return(data.frame(
      pathway = character(0),
      direction = character(0),
      n_pathway_genes = integer(0),
      n_sig_genes = integer(0),
      n_overlap = integer(0),
      expected_overlap = numeric(0),
      fold_enrichment = numeric(0),
      pvalue = numeric(0),
      genes = character(0),
      adj_pvalue = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # Sort by p-value
  results <- results[order(results$pvalue), ]
  rownames(results) <- NULL

  return(results)
}

#' Hypergeometric test for pathway enrichment
#' @keywords internal
.hypergeometricTest <- function(sig_genes, pathways, universe, min_genes, direction) {

  results <- lapply(names(pathways), function(pathway_name) {
    pathway_genes <- intersect(pathways[[pathway_name]], universe)

    if (length(pathway_genes) < min_genes) {
      return(NULL)
    }

    # Overlap
    overlap <- intersect(sig_genes, pathway_genes)

    if (length(overlap) == 0) {
      return(NULL)
    }

    # Hypergeometric test
    # x: number of white balls drawn (overlap)
    # m: number of white balls (pathway genes)
    # n: number of black balls (non-pathway genes)
    # k: number of balls drawn (significant genes)

    x <- length(overlap)
    m <- length(pathway_genes)
    n <- length(universe) - m
    k <- length(sig_genes)

    # P(X >= x)
    pval <- phyper(x - 1, m, n, k, lower.tail = FALSE)

    # Fold enrichment
    expected <- k * m / length(universe)
    fold_enrichment <- x / expected

    data.frame(
      pathway = pathway_name,
      direction = direction,
      n_pathway_genes = m,
      n_sig_genes = k,
      n_overlap = x,
      expected_overlap = round(expected, 2),
      fold_enrichment = fold_enrichment,
      pvalue = pval,
      genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results[!sapply(results, is.null)])

  # FIX: Check if results is empty or NULL
  if (is.null(results) || nrow(results) == 0) {
    # Return empty data.frame with proper structure
    return(data.frame(
      pathway = character(0),
      direction = character(0),
      n_pathway_genes = integer(0),
      n_sig_genes = integer(0),
      n_overlap = integer(0),
      expected_overlap = numeric(0),
      fold_enrichment = numeric(0),
      pvalue = numeric(0),
      genes = character(0),
      adj_pvalue = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  results$adj_pvalue <- p.adjust(results$pvalue, method = "BH")

  return(results)
}
