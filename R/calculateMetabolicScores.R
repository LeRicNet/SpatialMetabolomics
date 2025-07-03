#' Calculate metabolic pathway scores
#'
#' @param object SpatialMetabolic object
#' @param pathways List of pathway gene sets or pathway names
#' @param method Scoring method ("seurat", "ssgsea", "aucell", "mean")
#' @param scale Whether to scale scores (default: TRUE)
#' @param nbin Number of bins for Seurat scoring (default: 24)
#' @param ctrl Number of control genes for Seurat method (default: 100)
#' @param seed Random seed for reproducibility
#' @param BPPARAM BiocParallel parameters
#'
#' @return SpatialMetabolic object with metabolic scores
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom SingleCellExperiment logcounts logcounts<-
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom stats sd
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate scores for default pathways
#' spm <- calculateMetabolicScores(spm)
#'
#' # Calculate scores for specific pathways
#' spm <- calculateMetabolicScores(spm,
#'                                pathways = c("oxphos", "glycolysis"),
#'                                method = "seurat")
#'
#' # Use custom gene sets
#' custom_pathways <- list(
#'     MyPathway1 = c("Gene1", "Gene2", "Gene3"),
#'     MyPathway2 = c("Gene4", "Gene5", "Gene6")
#' )
#' spm <- calculateMetabolicScores(spm, pathways = custom_pathways)
#' }
setMethod("calculateMetabolicScores", "SpatialMetabolic",
          function(object,
                   pathways = c("oxphos", "glycolysis", "mito_biogenesis", "mito_dynamics"),
                   method = "seurat",
                   scale = TRUE,
                   nbin = 24,
                   ctrl = 100,
                   seed = 42,
                   verbose = TRUE,  # ADD THIS PARAMETER
                   BPPARAM = SerialParam()) {

            # Check if data is normalized
            if (all(logcounts(object) == 0)) {
              if (verbose) {
                message("No normalized data found. Running log normalization...")
              }
              object <- normalizeSpatial(object, verbose = verbose)
            }

            # Get pathway gene sets
            if (is.character(pathways)) {
              pathway_list <- .getDefaultPathways(pathways, object)
            } else if (is.list(pathways)) {
              pathway_list <- pathways
            } else {
              stop("pathways must be character vector or list of gene sets")
            }

            # Filter pathways to genes present in data
            pathway_list <- .filterPathwayGenes(pathway_list, object, verbose = verbose)

            if (length(pathway_list) == 0) {
              stop("No valid pathways after filtering for present genes")
            }

            # Store pathways
            metabolicPathways(object) <- pathway_list

            # Calculate scores based on method
            if (verbose) {
              message("Calculating metabolic scores using ", method, " method...")
            }

            set.seed(seed)

            if (method == "seurat") {
              scores <- .calculateSeuratScores(object, pathway_list, nbin, ctrl, BPPARAM)
            } else if (method == "mean") {
              scores <- .calculateMeanScores(object, pathway_list, BPPARAM)
            } else if (method == "ssgsea") {
              scores <- .calculateSSGSEAScores(object, pathway_list, BPPARAM)
            } else if (method == "aucell") {
              scores <- .calculateAUCellScores(object, pathway_list, BPPARAM)
            } else {
              stop("Method must be one of: seurat, mean, ssgsea, aucell")
            }

            # Scale if requested
            if (scale) {
              scores <- t(scale(t(scores)))
              # Replace any NAs from scaling with 0
              scores[is.na(scores)] <- 0
            }

            # Store scores
            metabolicScores(object) <- scores

            # Add to colData for easy access
            score_df <- as.data.frame(t(scores))
            colnames(score_df) <- paste0("score_", colnames(score_df))

            # Remove existing score columns if present
            existing_scores <- grep("^score_", colnames(colData(object)), value = TRUE)
            if (length(existing_scores) > 0) {
              colData(object) <- colData(object)[, !colnames(colData(object)) %in% existing_scores]
            }

            colData(object) <- cbind(colData(object), score_df)

            if (verbose) {
              message("Calculated scores for ", ncol(scores), " pathways")
            }

            return(object)
          })

#' Helper to get default pathways
#' @keywords internal
.getDefaultPathways <- function(pathway_names, object) {
  pathway_list <- list()

  # Determine species from gene names (simple heuristic)
  gene_names <- head(rownames(object), 100)
  if (mean(grepl("^[A-Z]", gene_names)) > 0.8) {
    species <- "human"
  } else {
    species <- "mouse"
  }

  for (pathway in pathway_names) {
    pathway <- tolower(pathway)
    if (pathway == "oxphos") {
      pathway_list[["OXPHOS"]] <- getOXPHOSGenes(species)
    } else if (pathway == "glycolysis") {
      pathway_list[["Glycolysis"]] <- getGlycolysisGenes(species)
    } else if (pathway == "mito_biogenesis") {
      pathway_list[["Mito_Biogenesis"]] <- getMitoBiogenesisGenes(species)
    } else if (pathway == "mito_dynamics") {
      pathway_list[["Mito_Dynamics"]] <- getMitoDynamicsGenes(species)
    } else if (pathway == "atp" || pathway == "atp_metabolism") {
      pathway_list[["ATP_Metabolism"]] <- getATPGenes(species)
    } else if (pathway == "svd" || pathway == "svd_signature") {
      pathway_list[["SVD_Signature"]] <- getSVDGenes(species)
    } else {
      warning("Unknown pathway: ", pathway)
    }
  }

  return(pathway_list)
}

#' Filter pathway genes to those present in data
#' @keywords internal
.filterPathwayGenes <- function(pathway_list, object, verbose = TRUE) {
  gene_universe <- rownames(object)

  filtered_list <- list()
  for (pathway_name in names(pathway_list)) {
    genes <- pathway_list[[pathway_name]]
    genes_present <- intersect(genes, gene_universe)

    if (length(genes_present) >= 3) {  # Minimum 3 genes per pathway
      filtered_list[[pathway_name]] <- genes_present
      coverage <- round(length(genes_present) / length(genes) * 100, 1)
      if (verbose) {
        message("  ", pathway_name, ": ", length(genes_present), "/",
                length(genes), " genes (", coverage, "% coverage)")
      }
    } else {
      if (verbose) {
        warning("Pathway ", pathway_name, " has too few genes (",
                length(genes_present), ") and will be skipped")
      }
    }
  }

  return(filtered_list)
}

#' Seurat-style scoring implementation
#' @keywords internal
#' @importFrom stats quantile
.calculateSeuratScores <- function(object, pathway_list, nbin, ctrl, BPPARAM) {
  expr_mat <- logcounts(object)

  # Calculate average expression for binning
  # Fix: Convert sparse matrix to regular matrix for rowMeans
  expr_mat <- as.matrix(expr_mat)
  avg_expr <- rowMeans(expr_mat)

  # Remove genes with no expression
  expressed_genes <- names(avg_expr)[avg_expr > 0]
  avg_expr <- avg_expr[expressed_genes]

  scores <- bplapply(pathway_list, function(genes) {
    genes_use <- intersect(genes, expressed_genes)

    if (length(genes_use) == 0) {
      return(rep(0, ncol(expr_mat)))
    }

    # Calculate pathway expression
    if (length(genes_use) == 1) {
      pathway_expr <- expr_mat[genes_use, ]
    } else {
      pathway_expr <- colMeans(expr_mat[genes_use, , drop = FALSE])
    }

    # Get control genes
    # Create bins based on average expression
    breaks <- quantile(avg_expr, probs = seq(0, 1, length.out = nbin + 1))
    breaks[1] <- breaks[1] - 1e-6  # Ensure lowest value is included
    expr_bins <- cut(avg_expr, breaks = breaks, include.lowest = TRUE)

    control_genes <- c()
    for (gene in genes_use) {
      if (gene %in% names(expr_bins)) {
        gene_bin <- expr_bins[gene]
        # Get all genes in the same bin
        bin_genes <- names(expr_bins)[expr_bins == gene_bin]
        # Remove pathway genes
        bin_genes <- setdiff(bin_genes, genes_use)

        if (length(bin_genes) >= ctrl) {
          # Sample control genes
          control_genes <- c(control_genes, sample(bin_genes, ctrl))
        } else if (length(bin_genes) > 0) {
          # Use all available genes if fewer than requested
          control_genes <- c(control_genes, bin_genes)
        }
      }
    }

    # Remove duplicates
    control_genes <- unique(control_genes)

    # Calculate control expression
    if (length(control_genes) > 1) {
      control_expr <- colMeans(expr_mat[control_genes, , drop = FALSE])
      score <- pathway_expr - control_expr
    } else {
      # If no good controls, just use pathway expression
      score <- pathway_expr
    }

    return(score)
  }, BPPARAM = BPPARAM)

  # Combine into matrix
  scores <- do.call(rbind, scores)
  rownames(scores) <- names(pathway_list)

  return(scores)
}

#' Simple mean expression scoring
#' @keywords internal
.calculateMeanScores <- function(object, pathway_list, BPPARAM) {
  expr_mat <- logcounts(object)

  scores <- bplapply(pathway_list, function(genes) {
    genes_use <- intersect(genes, rownames(expr_mat))

    if (length(genes_use) == 0) {
      return(rep(0, ncol(expr_mat)))
    }

    if (length(genes_use) == 1) {
      return(expr_mat[genes_use, ])
    } else {
      return(colMeans(expr_mat[genes_use, , drop = FALSE]))
    }
  }, BPPARAM = BPPARAM)

  scores <- do.call(rbind, scores)
  rownames(scores) <- names(pathway_list)

  return(scores)
}

#' Placeholder for ssGSEA scoring
#' @keywords internal
.calculateSSGSEAScores <- function(object, pathway_list, BPPARAM) {
  stop("ssGSEA scoring not yet implemented. Please use 'seurat' or 'mean' method.")
}

#' Placeholder for AUCell scoring
#' @keywords internal
.calculateAUCellScores <- function(object, pathway_list, BPPARAM) {
  stop("AUCell scoring not yet implemented. Please use 'seurat' or 'mean' method.")
}

#' Calculate pathway enrichment statistics
#'
#' @param object SpatialMetabolic object with calculated scores
#' @param group_by Column to group samples by
#' @return Data frame with enrichment statistics per group
#' @importFrom stats wilcox.test
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate enrichment by condition
#' enrichment <- calculatePathwayEnrichment(spm, group_by = "condition")
#' print(enrichment)
#' }
calculatePathwayEnrichment <- function(object, group_by = "condition") {

  if (ncol(metabolicScores(object)) == 0) {
    stop("No metabolic scores found. Run calculateMetabolicScores first.")
  }

  if (!group_by %in% colnames(colData(object))) {
    stop("Column '", group_by, "' not found in colData")
  }

  groups <- colData(object)[[group_by]]
  unique_groups <- unique(groups)

  if (length(unique_groups) < 2) {
    stop("Need at least 2 groups for enrichment analysis")
  }

  scores <- metabolicScores(object)
  results <- list()

  # Pairwise comparisons
  for (i in 1:(length(unique_groups)-1)) {
    for (j in (i+1):length(unique_groups)) {
      group1 <- unique_groups[i]
      group2 <- unique_groups[j]

      comparison <- paste0(group2, "_vs_", group1)

      pathway_results <- lapply(rownames(scores), function(pathway) {
        score1 <- scores[pathway, groups == group1]
        score2 <- scores[pathway, groups == group2]

        # Wilcoxon test
        if (length(score1) > 0 && length(score2) > 0 &&
            sd(c(score1, score2)) > 0) {
          wt <- wilcox.test(score2, score1)
          pval <- wt$p.value
        } else {
          pval <- 1
        }

        data.frame(
          pathway = pathway,
          comparison = comparison,
          mean_group1 = mean(score1),
          mean_group2 = mean(score2),
          log2FC = mean(score2) - mean(score1),
          pvalue = pval,
          stringsAsFactors = FALSE
        )
      })

      results[[comparison]] <- do.call(rbind, pathway_results)
      results[[comparison]]$adj_pvalue <- p.adjust(results[[comparison]]$pvalue,
                                                   method = "BH")
    }
  }

  # Combine all results
  combined_results <- do.call(rbind, results)
  rownames(combined_results) <- NULL

  return(combined_results)
}
