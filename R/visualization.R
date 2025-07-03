#' @title Visualization functions for SpatialMetabolic
#' @description Functions for creating publication-ready plots
#' @import ggplot2
#' @importFrom viridis scale_color_viridis_c scale_fill_viridis_c scale_fill_viridis_d
#' @importFrom scales percent

#' Plot spatial metabolic scores
#'
#' @param object SpatialMetabolic object
#' @param pathway Pathway name to plot
#' @param samples Sample IDs to plot (default: all)
#' @param ncol Number of columns for faceting
#' @param point_size Size of spatial points
#' @param color_scale Color scale ("viridis", "plasma", "inferno", "magma", "cividis")
#' @param limits Color scale limits (default: NULL for automatic)
#' @param na_color Color for NA values
#' @param legend_position Position of legend
#' @param title Plot title (NULL for automatic)
#'
#' @return ggplot object
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic plot
#' p <- plotSpatialMetabolicScore(spm, pathway = "OXPHOS")
#'
#' # Customize appearance
#' p <- plotSpatialMetabolicScore(spm,
#'                               pathway = "OXPHOS",
#'                               samples = c("WT1", "SVD1"),
#'                               color_scale = "plasma",
#'                               point_size = 2)
#' }
plotSpatialMetabolicScore <- function(object,
                                      pathway,
                                      samples = NULL,
                                      ncol = 2,
                                      point_size = 1.5,
                                      color_scale = "viridis",
                                      limits = NULL,
                                      na_color = "grey90",
                                      legend_position = "right",
                                      title = NULL) {

  # Check if scores exist
  if (nrow(metabolicScores(object)) == 0) {
    stop("No metabolic scores calculated. Run calculateMetabolicScores first.")
  }

  if (!pathway %in% rownames(metabolicScores(object))) {
    stop("Pathway '", pathway, "' not found in metabolic scores. Available pathways: ",
         paste(rownames(metabolicScores(object)), collapse = ", "))
  }

  # Get data
  plot_data <- data.frame(
    x = spatialCoords(object)[, 1],
    y = spatialCoords(object)[, 2],
    score = metabolicScores(object)[pathway, ],
    sample = colData(object)$sample_id,
    condition = comparisonGroups(object),
    stringsAsFactors = FALSE
  )

  # Filter samples if specified
  if (!is.null(samples)) {
    plot_data <- plot_data[plot_data$sample %in% samples, ]
    if (nrow(plot_data) == 0) {
      stop("No data found for specified samples")
    }
  }

  # Create facet labels
  plot_data$facet_label <- paste0(plot_data$sample, " (", plot_data$condition, ")")

  # Create plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = score)) +
    geom_point(size = point_size, shape = 16) +
    facet_wrap(~ facet_label, ncol = ncol) +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "grey80"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = legend_position,
      strip.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  # Add title
  if (is.null(title)) {
    p <- p + labs(title = paste(pathway, "Activity"), color = "Score")
  } else {
    p <- p + labs(title = title, color = "Score")
  }

  # Add color scale
  color_options <- c("viridis", "plasma", "inferno", "magma", "cividis")
  if (!color_scale %in% color_options) {
    warning("Unknown color scale. Using 'viridis'")
    color_scale <- "viridis"
  }

  p <- p + scale_color_viridis_c(option = color_scale,
                                 limits = limits,
                                 na.value = na_color)

  return(p)
}

#' Plot metabolic score comparison between conditions
#'
#' @param object SpatialMetabolic object
#' @param pathway Pathway name to plot
#' @param group_by Column to group by
#' @param plot_type Type of plot ("violin", "box", "density", "ridge")
#' @param add_points Whether to add individual points
#' @param point_alpha Alpha for points
#' @param colors Color palette for groups
#' @param add_stats Whether to add statistical annotations
#'
#' @return ggplot object
#' @importFrom stats wilcox.test
#' @export
plotMetabolicComparison <- function(object,
                                    pathway,
                                    group_by = "condition",
                                    plot_type = "violin",
                                    add_points = FALSE,
                                    point_alpha = 0.3,
                                    colors = NULL,
                                    add_stats = TRUE) {

  # Check inputs
  if (!pathway %in% rownames(metabolicScores(object))) {
    stop("Pathway '", pathway, "' not found")
  }

  if (!group_by %in% colnames(colData(object))) {
    stop("Column '", group_by, "' not found in colData")
  }

  # Get data
  plot_data <- data.frame(
    score = metabolicScores(object)[pathway, ],
    group = colData(object)[[group_by]],
    stringsAsFactors = FALSE
  )

  # Remove NAs
  plot_data <- plot_data[!is.na(plot_data$group), ]

  # Set colors
  n_groups <- length(unique(plot_data$group))
  if (is.null(colors)) {
    colors <- viridis::viridis(n_groups, option = "D")
  }

  # Create base plot
  p <- ggplot(plot_data, aes(x = group, y = score, fill = group))

  # Add plot type
  if (plot_type == "violin") {
    p <- p +
      geom_violin(alpha = 0.8, scale = "width") +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8,
                   outlier.shape = NA)
  } else if (plot_type == "box") {
    p <- p + geom_boxplot(alpha = 0.8, outlier.shape = NA)
  } else if (plot_type == "density") {
    p <- ggplot(plot_data, aes(x = score, fill = group)) +
      geom_density(alpha = 0.7) +
      facet_wrap(~ group, ncol = 1, scales = "free_y")
  } else if (plot_type == "ridge") {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("Package 'ggridges' needed for ridge plots")
    }
    p <- ggplot(plot_data, aes(x = score, y = group, fill = group)) +
      ggridges::geom_density_ridges(alpha = 0.8)
  }

  # Add points if requested
  if (add_points && plot_type %in% c("violin", "box")) {
    p <- p +
      geom_jitter(width = 0.2, alpha = point_alpha, size = 0.5,
                  shape = 16, color = "black")
  }

  # Styling
  p <- p +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = paste(pathway, "Activity by", group_by),
      y = "Metabolic Score",
      x = ""
    ) +
    scale_fill_manual(values = colors)

  # Add statistics if requested and only 2 groups
  if (add_stats && plot_type %in% c("violin", "box")) {
    groups <- unique(plot_data$group)
    if (length(groups) == 2) {
      # Wilcoxon test
      test_result <- wilcox.test(
        plot_data$score[plot_data$group == groups[1]],
        plot_data$score[plot_data$group == groups[2]]
      )

      # Format p-value
      if (test_result$p.value < 0.001) {
        p_text <- "p < 0.001"
      } else {
        p_text <- paste("p =", format(test_result$p.value, digits = 3))
      }

      # Add annotation
      y_max <- max(plot_data$score) * 1.1
      p <- p +
        annotate("segment", x = 1, xend = 2, y = y_max, yend = y_max) +
        annotate("text", x = 1.5, y = y_max * 1.02, label = p_text,
                 size = 3.5)
    }
  }

  return(p)
}

#' Plot spatial gradients
#'
#' @param object SpatialMetabolic object
#' @param top_n Number of top features to plot
#' @param result_type Which gradient results to use ("scores" or "genes")
#' @param show_names Whether to show feature names
#' @param name_size Size of feature names
#'
#' @return ggplot object
#' @export
plotSpatialGradients <- function(object,
                                 top_n = 10,
                                 result_type = "scores",
                                 show_names = TRUE,
                                 name_size = 3) {

  # Get gradient results
  result_name <- paste0("spatial_gradients_", result_type)

  if (!result_name %in% names(analysisResults(object))) {
    stop("No spatial gradient results found for ", result_type,
         ". Run detectMetabolicGradients first.")
  }

  grad_results <- analysisResults(object)[[result_name]]

  # Filter to valid results
  grad_results <- grad_results[!is.na(grad_results$adj_pvalue), ]

  # Get top features
  top_features <- head(grad_results[order(grad_results$adj_pvalue), ], top_n)

  # Prepare data for plotting
  top_features$neg_log10_p <- -log10(top_features$adj_pvalue)
  top_features$significant <- top_features$adj_pvalue < 0.05

  # Reorder by significance
  top_features$feature <- factor(top_features$feature,
                                 levels = rev(top_features$feature))

  # Create plot
  p <- ggplot(top_features, aes(x = neg_log10_p, y = feature, fill = significant)) +
    geom_col() +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"),
                      labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = ifelse(show_names, 9, 0)),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "top",
      legend.title = element_blank()
    ) +
    labs(
      title = paste("Top", top_n, "Spatial Metabolic Patterns"),
      x = "-log10(Adjusted P-value)",
      y = ""
    )

  # Add Moran's I values if available
  if ("morans_i" %in% colnames(top_features)) {
    p <- p +
      geom_text(aes(label = sprintf("I=%.2f", morans_i)),
                hjust = -0.1, size = name_size)
  }

  return(p)
}

#' Plot metabolic heatmap
#'
#' @param object SpatialMetabolic object
#' @param features Features to include in heatmap (NULL for metabolic scores)
#' @param group_by Column to group samples by
#' @param scale_rows Whether to scale rows
#' @param show_row_names Whether to show row names
#' @param show_column_names Whether to show column names
#' @param cluster_rows Whether to cluster rows
#' @param cluster_cols Whether to cluster columns
#' @param annotation_colors List of colors for annotations
#' @param color_palette Color palette for heatmap
#' @param ... Additional arguments passed to pheatmap
#'
#' @return pheatmap object or ggplot object
#' @importFrom stats hclust dist
#' @export
plotMetabolicHeatmap <- function(object,
                                 features = NULL,
                                 group_by = c("condition", "sample_id"),
                                 scale_rows = TRUE,
                                 show_row_names = TRUE,
                                 show_column_names = FALSE,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 annotation_colors = NULL,
                                 color_palette = "RdBu",
                                 ...) {

  # Check if pheatmap is available
  use_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

  # Get features
  if (is.null(features)) {
    if (nrow(metabolicScores(object)) > 0) {
      mat <- metabolicScores(object)
      # Calculate per-sample means
      samples <- unique(colData(object)$sample_id)
      sample_means <- sapply(samples, function(s) {
        idx <- colData(object)$sample_id == s
        rowMeans(mat[, idx, drop = FALSE])
      })
      mat <- sample_means
    } else {
      stop("No features specified and no metabolic scores available")
    }
  } else {
    if (is.character(features)) {
      mat <- logcounts(object)[features, , drop = FALSE]
    } else {
      stop("features must be character vector or NULL")
    }
  }

  # Scale if requested
  if (scale_rows) {
    mat <- t(scale(t(mat)))
    # Handle any NAs from scaling
    mat[is.na(mat)] <- 0
  }

  # Get annotation
  if (ncol(mat) == length(unique(colData(object)$sample_id))) {
    # Per-sample data
    sample_info <- unique(colData(object)[, c("sample_id", group_by), drop = FALSE])
    # Fix duplicate row names issue
    sample_info <- sample_info[!duplicated(sample_info$sample_id), ]
    rownames(sample_info) <- sample_info$sample_id
    col_anno <- sample_info[colnames(mat), group_by, drop = FALSE]
  } else {
    # Per-spot data
    col_anno <- colData(object)[, group_by, drop = FALSE]
  }

  # Convert col_anno to data.frame if it isn't already
  if (!is.data.frame(col_anno)) {
    col_anno <- as.data.frame(col_anno)
    colnames(col_anno) <- group_by
  }

  if (use_pheatmap) {
    # Use pheatmap

    # Set up colors
    if (is.null(annotation_colors)) {
      annotation_colors <- list()
      for (col in group_by) {
        if (col %in% colnames(col_anno)) {
          n_levels <- length(unique(col_anno[[col]]))
          annotation_colors[[col]] <- viridis::viridis(n_levels, option = "D")
          names(annotation_colors[[col]]) <- unique(col_anno[[col]])
        }
      }
    }

    # Color palette
    if (color_palette == "RdBu") {
      colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100)
    } else if (color_palette == "viridis") {
      colors <- viridis::viridis(100)
    } else {
      colors <- colorRampPalette(c("blue", "white", "red"))(100)
    }

    # Create heatmap
    p <- pheatmap::pheatmap(
      mat,
      annotation_col = col_anno,
      annotation_colors = annotation_colors,
      scale = "none",  # Already scaled if requested
      show_rownames = show_row_names,
      show_colnames = show_column_names,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      color = colors,
      border_color = NA,
      ...
    )

    return(p)

  } else {
    # Fallback to ggplot2 implementation
    message("pheatmap not available, using ggplot2 implementation")

    # Prepare data for ggplot
    mat_long <- reshape2::melt(as.matrix(mat))
    colnames(mat_long) <- c("Feature", "Sample", "Value")

    # Add annotations
    mat_long <- merge(mat_long, col_anno, by.x = "Sample", by.y = "row.names")

    # Create plot
    p <- ggplot(mat_long, aes(x = Sample, y = Feature, fill = Value)) +
      geom_tile() +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = ifelse(show_column_names, 8, 0)),
        axis.text.y = element_text(size = ifelse(show_row_names, 8, 0)),
        axis.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank()
      )

    # Add color scale
    if (scale_rows) {
      p <- p + scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                    midpoint = 0, name = "Z-score")
    } else {
      p <- p + scale_fill_viridis_c(name = "Expression")
    }

    # Add faceting for groups if only one group variable
    if (length(group_by) == 1) {
      p <- p + facet_grid(. ~ get(group_by[1]), scales = "free_x", space = "free")
    }

    return(p)
  }
}

#' Create spatial feature plot for any gene or score
#'
#' @param object SpatialMetabolic object
#' @param features Character vector of features to plot
#' @param ncol Number of columns for multiple features
#' @param samples Samples to include
#' @param point_size Size of points
#' @param color_scale Color scale to use
#' @param sync_scales Whether to synchronize color scales across features
#'
#' @return ggplot object
#' @export
plotSpatialFeatures <- function(object,
                                features,
                                ncol = 2,
                                samples = NULL,
                                point_size = 1.5,
                                color_scale = "viridis",
                                sync_scales = TRUE) {

  # Determine feature types
  gene_features <- features[features %in% rownames(object)]
  score_features <- features[features %in% rownames(metabolicScores(object))]

  missing <- setdiff(features, c(gene_features, score_features))
  if (length(missing) > 0) {
    warning("Features not found: ", paste(missing, collapse = ", "))
  }

  # Collect all data
  plot_list <- list()

  # Gene features
  if (length(gene_features) > 0) {
    gene_data <- lapply(gene_features, function(gene) {
      data.frame(
        x = spatialCoords(object)[, 1],
        y = spatialCoords(object)[, 2],
        value = logcounts(object)[gene, ],
        feature = gene,
        feature_type = "Gene",
        sample = colData(object)$sample_id,
        condition = comparisonGroups(object),
        stringsAsFactors = FALSE
      )
    })
    plot_list <- c(plot_list, gene_data)
  }

  # Score features
  if (length(score_features) > 0) {
    score_data <- lapply(score_features, function(score) {
      data.frame(
        x = spatialCoords(object)[, 1],
        y = spatialCoords(object)[, 2],
        value = metabolicScores(object)[score, ],
        feature = score,
        feature_type = "Pathway",
        sample = colData(object)$sample_id,
        condition = comparisonGroups(object),
        stringsAsFactors = FALSE
      )
    })
    plot_list <- c(plot_list, score_data)
  }

  # Combine data
  plot_data <- do.call(rbind, plot_list)

  # Filter samples
  if (!is.null(samples)) {
    plot_data <- plot_data[plot_data$sample %in% samples, ]
  }

  # Calculate limits if syncing scales
  if (sync_scales) {
    value_range <- range(plot_data$value, na.rm = TRUE)
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = value)) +
    geom_point(size = point_size, shape = 16) +
    facet_wrap(~ feature + feature_type, ncol = ncol, scales = "free") +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right"
    )

  # Add color scale
  if (sync_scales) {
    p <- p + scale_color_viridis_c(option = color_scale, limits = value_range)
  } else {
    p <- p + scale_color_viridis_c(option = color_scale)
  }

  return(p)
}

#' Plot metabolic trajectory analysis
#'
#' @param object SpatialMetabolic object
#' @param trajectory_scores Trajectory scores from trajectory analysis
#' @param color_by Variable to color by ("pseudotime", "condition", "pathway")
#' @param facet_by Variable to facet by
#' @param smooth_method Smoothing method for trajectory
#'
#' @return ggplot object
#' @export
plotMetabolicTrajectory <- function(object,
                                    trajectory_scores,
                                    color_by = "pseudotime",
                                    facet_by = NULL,
                                    smooth_method = "loess") {

  # This is a placeholder for trajectory analysis visualization
  # Would require trajectory analysis implementation first

  stop("Metabolic trajectory analysis not yet implemented")
}

#' Create summary plot of all metabolic pathways
#'
#' @param object SpatialMetabolic object
#' @param group_by Grouping variable
#' @param plot_type Type of summary plot ("box", "heatmap", "radar")
#' @param show_significance Whether to show significance annotations
#'
#' @return ggplot object
#' @export
plotMetabolicSummary <- function(object,
                                 group_by = "condition",
                                 plot_type = "box",
                                 show_significance = TRUE) {

  if (nrow(metabolicScores(object)) == 0) {
    stop("No metabolic scores found")
  }

  # Prepare data
  scores_long <- reshape2::melt(metabolicScores(object))
  colnames(scores_long) <- c("Pathway", "Cell", "Score")

  # Add metadata
  scores_long$Group <- rep(colData(object)[[group_by]],
                           each = nrow(metabolicScores(object)))

  if (plot_type == "box") {
    p <- ggplot(scores_long, aes(x = Pathway, y = Score, fill = Group)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      ) +
      scale_fill_viridis_d() +
      labs(title = "Metabolic Pathway Activity Summary")

    # Add significance if requested
    if (show_significance && length(unique(scores_long$Group)) == 2) {
      # Calculate p-values for each pathway
      pathways <- unique(scores_long$Pathway)
      groups <- unique(scores_long$Group)

      for (i in seq_along(pathways)) {  # FIXED: Proper R syntax
        pathway <- pathways[i]
        pathway_data <- scores_long[scores_long$Pathway == pathway, ]

        p_val <- wilcox.test(
          pathway_data$Score[pathway_data$Group == groups[1]],
          pathway_data$Score[pathway_data$Group == groups[2]]
        )$p.value

        if (p_val < 0.05) {
          y_pos <- max(pathway_data$Score) * 1.1

          if (p_val < 0.001) {
            label <- "***"
          } else if (p_val < 0.01) {
            label <- "**"
          } else {
            label <- "*"
          }

          p <- p + annotate("text", x = i, y = y_pos, label = label)
        }
      }
    }

  } else if (plot_type == "heatmap") {
    # Calculate mean scores per group and pathway
    mean_scores <- aggregate(Score ~ Pathway + Group, scores_long, mean)
    mean_matrix <- reshape2::dcast(mean_scores, Pathway ~ Group, value.var = "Score")
    rownames(mean_matrix) <- mean_matrix$Pathway
    mean_matrix$Pathway <- NULL

    # Convert to long format for ggplot
    heatmap_data <- reshape2::melt(as.matrix(mean_matrix))
    colnames(heatmap_data) <- c("Pathway", "Group", "Score")

    p <- ggplot(heatmap_data, aes(x = Group, y = Pathway, fill = Score)) +
      geom_tile() +
      scale_fill_viridis_c() +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()
      ) +
      labs(title = "Mean Metabolic Pathway Activity")

  } else if (plot_type == "radar") {
    stop("Radar plots not yet implemented")
  }

  return(p)
}

#' Volcano plot for differential expression results
#'
#' @param de_results Results from compareMetabolicStates
#' @param pval_cutoff P-value cutoff
#' @param fc_cutoff Fold change cutoff (absolute value)
#' @param label_top Number of top features to label
#' @param colors Colors for down, non-sig, and up
#'
#' @return ggplot object
#' @importFrom ggrepel geom_text_repel
#' @export
plotVolcano <- function(de_results,
                        pval_cutoff = 0.05,
                        fc_cutoff = 0.5,
                        label_top = 10,
                        colors = c("blue", "grey", "red")) {

  # Prepare data
  de_results$neg_log10_p <- -log10(de_results$adj_pvalue)
  de_results$neg_log10_p[is.infinite(de_results$neg_log10_p)] <-
    max(de_results$neg_log10_p[!is.infinite(de_results$neg_log10_p)]) * 1.1

  # Classify points
  de_results$significance <- "Not significant"
  de_results$significance[de_results$adj_pvalue < pval_cutoff &
                            de_results$log2FC > fc_cutoff] <- "Up"
  de_results$significance[de_results$adj_pvalue < pval_cutoff &
                            de_results$log2FC < -fc_cutoff] <- "Down"
  de_results$significance <- factor(de_results$significance,
                                    levels = c("Down", "Not significant", "Up"))

  # Create plot
  p <- ggplot(de_results, aes(x = log2FC, y = neg_log10_p, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(pval_cutoff),
               linetype = "dashed", color = "grey50") +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    ) +
    labs(
      title = "Differential Expression Volcano Plot",
      x = "log2 Fold Change",
      y = "-log10(Adjusted P-value)"
    )

  # Add labels for top features
  if (label_top > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    # Get top up and down features
    top_up <- de_results[de_results$significance == "Up", ]
    top_up <- top_up[order(top_up$adj_pvalue), ]
    top_up <- head(top_up, ceiling(label_top / 2))

    top_down <- de_results[de_results$significance == "Down", ]
    top_down <- top_down[order(top_down$adj_pvalue), ]
    top_down <- head(top_down, floor(label_top / 2))

    top_features <- rbind(top_up, top_down)

    if (nrow(top_features) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = top_features,
        aes(label = feature),
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50",
        max.overlaps = 20
      )
    }
  }

  return(p)
}
