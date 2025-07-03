#' Additional visualization functions for SpatialMetabolics
#' @import ggplot2
#' @importFrom cowplot theme_cowplot

#' Create spatial expression heatmap grid
#'
#' @param object SpatialMetabolic object
#' @param features Features to plot in grid
#' @param sample Sample to plot
#' @param ncol Number of columns
#' @param point_size Point size
#' @param sync_scale Synchronize color scales
#' @param color_palette Color palette
#'
#' @return ggplot object
#' @export
plotSpatialGrid <- function(object,
                            features,
                            sample = NULL,
                            ncol = 4,
                            point_size = 0.8,
                            sync_scale = TRUE,
                            color_palette = "viridis") {

  # Get sample
  if (is.null(sample)) {
    sample <- unique(colData(object)$sample_id)[1]
  }

  # Subset to sample
  sample_idx <- colData(object)$sample_id == sample

  # Collect data for all features
  plot_data_list <- lapply(features, function(f) {
    # Check if gene or score
    if (f %in% rownames(object)) {
      values <- logcounts(object)[f, sample_idx]
      feature_type <- "Gene"
    } else if (f %in% rownames(metabolicScores(object))) {
      values <- metabolicScores(object)[f, sample_idx]
      feature_type <- "Score"
    } else {
      return(NULL)
    }

    data.frame(
      x = spatialCoords(object)[sample_idx, 1],
      y = spatialCoords(object)[sample_idx, 2],
      value = values,
      feature = f,
      feature_type = feature_type,
      stringsAsFactors = FALSE
    )
  })

  # Remove NULL entries
  plot_data_list <- plot_data_list[!sapply(plot_data_list, is.null)]

  if (length(plot_data_list) == 0) {
    stop("No valid features found")
  }

  # Combine data
  plot_data <- do.call(rbind, plot_data_list)

  # Calculate shared scale if requested
  if (sync_scale) {
    value_limits <- range(plot_data$value, na.rm = TRUE)
  } else {
    value_limits <- NULL
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = value)) +
    geom_point(size = point_size) +
    facet_wrap(~ feature, ncol = ncol) +
    coord_fixed() +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.spacing = unit(0.5, "lines"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1, "cm")
    ) +
    labs(title = paste("Spatial Expression:", sample))

  # Add color scale
  if (color_palette == "viridis") {
    p <- p + scale_color_viridis_c(limits = value_limits, name = "Expression")
  } else if (color_palette == "plasma") {
    p <- p + scale_color_viridis_c(option = "plasma", limits = value_limits, name = "Expression")
  } else if (color_palette == "RdBu") {
    p <- p + scale_color_gradient2(low = "blue", mid = "white", high = "red",
                                   limits = value_limits, name = "Expression")
  }

  return(p)
}

#' Plot metabolic correlation network
#'
#' @param object SpatialMetabolic object
#' @param min_correlation Minimum correlation to show edge
#' @param layout Network layout algorithm
#' @param node_size Node size
#' @param edge_width_range Range for edge widths
#' @param colors Node colors
#'
#' @return ggplot object
#' @export
plotMetabolicNetwork <- function(object,
                                 min_correlation = 0.3,
                                 layout = "fr",
                                 node_size = 5,
                                 edge_width_range = c(0.5, 3),
                                 colors = NULL) {

  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    stop("Packages 'igraph' and 'ggraph' required for network plots")
  }

  # Calculate correlations between pathways
  if (nrow(metabolicScores(object)) < 2) {
    stop("Need at least 2 pathways for network")
  }

  cor_mat <- cor(t(metabolicScores(object)))

  # Create network
  adj_mat <- abs(cor_mat)
  adj_mat[adj_mat < min_correlation] <- 0
  diag(adj_mat) <- 0

  g <- igraph::graph_from_adjacency_matrix(
    adj_mat,
    mode = "undirected",
    weighted = TRUE
  )

  if (igraph::ecount(g) == 0) {
    stop("No edges above correlation threshold")
  }

  # Create layout
  if (layout == "fr") {
    lay <- igraph::layout_with_fr(g)
  } else if (layout == "circle") {
    lay <- igraph::layout_in_circle(g)
  } else {
    lay <- igraph::layout_nicely(g)
  }

  # Prepare for ggraph
  g_df <- ggraph::create_layout(g, layout = lay)

  # Add node attributes
  g_df$pathway <- rownames(cor_mat)

  # Create plot
  p <- ggraph::ggraph(g_df) +
    ggraph::geom_edge_link(aes(width = weight),
                           alpha = 0.6,
                           color = "grey50") +
    ggraph::geom_node_point(size = node_size,
                            aes(color = pathway)) +
    ggraph::geom_node_text(aes(label = pathway),
                           repel = TRUE,
                           size = 3) +
    ggraph::scale_edge_width(range = edge_width_range) +
    theme_void() +
    theme(legend.position = "none") +
    labs(title = "Metabolic Pathway Correlation Network")

  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  } else {
    p <- p + scale_color_viridis_d()
  }

  return(p)
}

#' Plot sample-to-sample correlation heatmap
#'
#' @param object SpatialMetabolic object
#' @param use_pathways Whether to use pathway scores or genes
#' @param method Correlation method
#' @param cluster Whether to cluster samples
#'
#' @return pheatmap or ggplot object
#' @export
plotSampleCorrelation <- function(object,
                                  use_pathways = TRUE,
                                  method = "pearson",
                                  cluster = TRUE) {

  # Get data matrix
  if (use_pathways && nrow(metabolicScores(object)) > 0) {
    # Average by sample
    samples <- unique(colData(object)$sample_id)
    data_mat <- sapply(samples, function(s) {
      idx <- colData(object)$sample_id == s
      rowMeans(metabolicScores(object)[, idx, drop = FALSE])
    })
  } else {
    # Use top variable genes
    var_genes <- rownames(object)[rowData(object)$variable_feature]
    if (length(var_genes) == 0) {
      var_genes <- head(rownames(object), 500)
    }

    samples <- unique(colData(object)$sample_id)
    data_mat <- sapply(samples, function(s) {
      idx <- colData(object)$sample_id == s
      rowMeans(logcounts(object)[var_genes, idx, drop = FALSE])
    })
  }

  # Calculate correlation
  cor_mat <- cor(data_mat, method = method)

  # Get sample annotations
  sample_meta <- unique(colData(object)[, c("sample_id", "condition")])
  # Fix duplicate row names issue
  sample_meta <- sample_meta[!duplicated(sample_meta$sample_id), ]
  rownames(sample_meta) <- sample_meta$sample_id
  sample_meta <- sample_meta[colnames(cor_mat), , drop = FALSE]

  # Ensure it's a proper data.frame
  if (!is.data.frame(sample_meta)) {
    sample_meta <- data.frame(
      condition = sample_meta,
      row.names = colnames(cor_mat)
    )
  }

  # Create heatmap
  if (requireNamespace("pheatmap", quietly = TRUE)) {

    # Annotation colors
    ann_colors <- list(
      condition = c(WT = "blue", SVD = "red", Disease = "red")
    )

    # Only use pheatmap if we have valid annotations
    if (nrow(sample_meta) > 0 && ncol(sample_meta) > 0) {
      p <- pheatmap::pheatmap(
        cor_mat,
        annotation_col = sample_meta,
        annotation_colors = ann_colors,
        cluster_rows = cluster,
        cluster_cols = cluster,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-1, 1, length.out = 101),
        main = "Sample Correlation",
        display_numbers = TRUE,
        number_format = "%.2f"
      )
    } else {
      # Fallback without annotations
      p <- pheatmap::pheatmap(
        cor_mat,
        cluster_rows = cluster,
        cluster_cols = cluster,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-1, 1, length.out = 101),
        main = "Sample Correlation",
        display_numbers = TRUE,
        number_format = "%.2f"
      )
    }
  }

  return(p)
}

#' Create comprehensive QC plot
#'
#' @param object SpatialMetabolic object
#' @param metrics QC metrics to plot
#' @param ncol Number of columns
#'
#' @return ggplot object
#' @export
plotQCSpatial <- function(object,
                          metrics = c("nCount_RNA", "nFeature_RNA",
                                      "percent_mt", "percent_ribo"),
                          ncol = 2) {

  # Check which metrics are available
  available_metrics <- intersect(metrics, colnames(colData(object)))

  if (length(available_metrics) == 0) {
    stop("No QC metrics found. Run calculateQCMetrics first.")
  }

  # Prepare data
  plot_data_list <- lapply(available_metrics, function(metric) {
    data.frame(
      x = spatialCoords(object)[, 1],
      y = spatialCoords(object)[, 2],
      value = colData(object)[[metric]],
      metric = metric,
      sample = colData(object)$sample_id,
      stringsAsFactors = FALSE
    )
  })

  plot_data <- do.call(rbind, plot_data_list)

  # Create plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    facet_grid(metric ~ sample, scales = "free") +
    coord_fixed() +
    scale_color_viridis_c() +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(size = 10)
    ) +
    labs(title = "Spatial QC Metrics")

  return(p)
}

#' Plot metabolic state transitions
#'
#' @param object SpatialMetabolic object
#' @param pathway1 First pathway
#' @param pathway2 Second pathway
#' @param group_by Grouping variable
#' @param add_ellipse Add confidence ellipses
#' @param point_alpha Point transparency
#'
#' @return ggplot object
#' @export
plotMetabolicScatter <- function(object,
                                 pathway1,
                                 pathway2,
                                 group_by = "condition",
                                 add_ellipse = TRUE,
                                 point_alpha = 0.6) {

  # Check pathways exist
  if (!all(c(pathway1, pathway2) %in% rownames(metabolicScores(object)))) {
    stop("Pathways not found in metabolic scores")
  }

  # Prepare data
  plot_data <- data.frame(
    x = metabolicScores(object)[pathway1, ],
    y = metabolicScores(object)[pathway2, ],
    group = colData(object)[[group_by]],
    stringsAsFactors = FALSE
  )

  # Create plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_point(alpha = point_alpha, size = 1) +
    theme_minimal() +
    labs(
      title = paste("Metabolic State:", pathway1, "vs", pathway2),
      x = paste(pathway1, "Score"),
      y = paste(pathway2, "Score")
    ) +
    scale_color_viridis_d()

  # Add ellipses if requested
  if (add_ellipse && requireNamespace("ggforce", quietly = TRUE)) {
    p <- p + ggforce::geom_mark_ellipse(
      aes(fill = group),
      alpha = 0.1,
      expand = unit(0.5, "mm")
    )
  } else if (add_ellipse) {
    p <- p + stat_ellipse(type = "norm", level = 0.95)
  }

  return(p)
}

#' Create animated spatial plot
#'
#' @param object SpatialMetabolic object
#' @param pathways Pathways to animate through
#' @param sample Sample to plot
#' @param fps Frames per second
#' @param duration Duration in seconds
#'
#' @return gganimate object
#' @export
animateSpatialMetabolics <- function(object,
                                     pathways = NULL,
                                     sample = NULL,
                                     fps = 10,
                                     duration = 10) {

  if (!requireNamespace("gganimate", quietly = TRUE)) {
    stop("Package 'gganimate' required for animations")
  }

  # Get pathways
  if (is.null(pathways)) {
    pathways <- rownames(metabolicScores(object))
  }

  # Get sample
  if (is.null(sample)) {
    sample <- unique(colData(object)$sample_id)[1]
  }

  # Subset to sample
  sample_idx <- colData(object)$sample_id == sample

  # Prepare data
  plot_data_list <- lapply(pathways, function(p) {
    data.frame(
      x = spatialCoords(object)[sample_idx, 1],
      y = spatialCoords(object)[sample_idx, 2],
      value = metabolicScores(object)[p, sample_idx],
      pathway = p,
      stringsAsFactors = FALSE
    )
  })

  plot_data <- do.call(rbind, plot_data_list)

  # Create animated plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = value)) +
    geom_point(size = 2) +
    scale_color_viridis_c() +
    coord_fixed() +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = '{closest_state}', color = "Score") +
    gganimate::transition_states(pathway,
                                 transition_length = 2,
                                 state_length = 1) +
    gganimate::enter_fade() +
    gganimate::exit_fade()

  # Animate
  anim <- gganimate::animate(p,
                             fps = fps,
                             duration = duration,
                             width = 400,
                             height = 400)

  return(anim)
}

#' Create summary report plot
#'
#' @param object SpatialMetabolic object
#' @param save_to File path to save report
#'
#' @return Combined plot object
#' @export
createMetabolicReport <- function(object, save_to = NULL) {

  # Check data availability
  if (nrow(metabolicScores(object)) == 0) {
    stop("No metabolic scores calculated")
  }

  # 1. Sample correlation
  p1 <- plotSampleCorrelation(object, cluster = FALSE)

  # 2. Top pathways by condition
  pathway_means <- aggregate(
    t(metabolicScores(object)),
    by = list(condition = comparisonGroups(object)),
    FUN = mean
  )

  pathway_long <- reshape2::melt(pathway_means, id.vars = "condition")

  p2 <- ggplot(pathway_long, aes(x = variable, y = value, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mean Pathway Activity",
         x = "", y = "Score") +
    scale_fill_viridis_d()

  # 3. PCA if available
  if ("PCA" %in% reducedDimNames(object)) {
    pca_data <- data.frame(
      PC1 = reducedDim(object, "PCA")[, 1],
      PC2 = reducedDim(object, "PCA")[, 2],
      condition = comparisonGroups(object),
      sample = colData(object)$sample_id
    )

    p3 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(alpha = 0.6) +
      stat_ellipse() +
      theme_minimal() +
      scale_color_viridis_d() +
      labs(title = "PCA: Metabolic Space")
  } else {
    p3 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "PCA not calculated") +
      theme_void()
  }

  # 4. Summary statistics
  n_samples <- length(unique(colData(object)$sample_id))
  n_conditions <- length(unique(comparisonGroups(object)))
  n_pathways <- nrow(metabolicScores(object))
  n_spots <- ncol(object)

  summary_text <- paste(
    "Dataset Summary:",
    paste("- Samples:", n_samples),
    paste("- Conditions:", n_conditions),
    paste("- Spots:", format(n_spots, big.mark = ",")),
    paste("- Pathways analyzed:", n_pathways),
    sep = "\n"
  )

  p4 <- ggplot() +
    annotate("text", x = 0.1, y = 0.5, label = summary_text,
             hjust = 0, vjust = 0.5, size = 4) +
    theme_void() +
    labs(title = "Analysis Summary")

  # Combine
  report <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "SpatialMetabolics Analysis Report",
      subtitle = paste("Generated:", Sys.Date())
    )

  # Save if requested
  if (!is.null(save_to)) {
    ggsave(save_to, report, width = 12, height = 10)
    message("Report saved to: ", save_to)
  }

  return(report)
}
