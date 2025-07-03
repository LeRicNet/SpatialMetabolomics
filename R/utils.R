#' @title Utility functions for SpatialMetabolic
#' @description Helper functions for data processing and analysis

#' Normalize spatial data
#'
#' @param object SpatialMetabolic object
#' @param method Normalization method ("logNormalize", "SCTransform", "TMM", "scran")
#' @param scale_factor Scale factor for log normalization
#' @param verbose Print messages
#'
#' @return Normalized SpatialMetabolic object
#' @importFrom SingleCellExperiment logcounts logcounts<- counts
#' @importFrom Matrix colSums rowSums
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic log normalization
#' spm <- normalizeSpatial(spm)
#'
#' # Use different scale factor
#' spm <- normalizeSpatial(spm, scale_factor = 1e6)
#' }
normalizeSpatial <- function(object,
                             method = "logNormalize",
                             scale_factor = 10000,
                             verbose = TRUE) {

  if (verbose) {
    message("Normalizing data using ", method, " method...")
  }

  if (method == "logNormalize") {
    # Get counts
    count_mat <- counts(object)

    # Calculate library sizes
    lib_sizes <- colSums(count_mat)

    # Avoid division by zero
    lib_sizes[lib_sizes == 0] <- 1

    # Normalize
    norm_counts <- t(t(count_mat) / lib_sizes * scale_factor)

    # Log transform
    log_counts <- log1p(norm_counts)

    # Store
    logcounts(object) <- log_counts

  } else if (method == "SCTransform") {
    # Requires Seurat integration
    stop("SCTransform normalization requires Seurat. ",
         "Convert to Seurat object first or use 'logNormalize'")

  } else if (method == "TMM") {
    # TMM normalization (requires edgeR)
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package 'edgeR' required for TMM normalization")
    }

    # Create DGEList
    dge <- edgeR::DGEList(counts = counts(object))

    # Calculate normalization factors
    dge <- edgeR::calcNormFactors(dge, method = "TMM")

    # Get normalized counts
    norm_counts <- edgeR::cpm(dge, log = FALSE)
    log_counts <- log1p(norm_counts)

    logcounts(object) <- log_counts

  } else if (method == "scran") {
    # scran normalization
    if (!requireNamespace("scran", quietly = TRUE)) {
      stop("Package 'scran' required for scran normalization")
    }

    # Calculate size factors
    size_factors <- scran::calculateSumFactors(object)

    # Normalize
    object <- scran::logNormCounts(object, size.factors = size_factors)

  } else {
    stop("Method must be one of: logNormalize, SCTransform, TMM, scran")
  }

  if (verbose) {
    message("Normalization complete")
  }

  return(object)
}

#' Preprocess SpatialMetabolic object
#'
#' @param object SpatialMetabolic object
#' @param normalize_method Normalization method
#' @param n_variable_features Number of variable features
#' @param n_pcs Number of principal components
#' @param verbose Print messages
#'
#' @return Preprocessed object
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
#' @importFrom SummarizedExperiment rowData rowData<-
#' @export
preprocessSpatial <- function(object,
                              normalize_method = "logNormalize",
                              n_variable_features = 2000,
                              n_pcs = 30,
                              verbose = TRUE) {

  if (verbose) {
    message("Preprocessing spatial data...")
  }

  # Normalize if not already done
  if (all(logcounts(object) == 0)) {
    object <- normalizeSpatial(object, method = normalize_method, verbose = verbose)
  }

  # Find variable features
  if (n_variable_features > 0) {
    object <- findVariableFeatures(object,
                                   n_features = n_variable_features,
                                   verbose = verbose)
  }

  # Run PCA
  if (n_pcs > 0) {
    object <- runPCA(object, n_components = n_pcs, verbose = verbose)
  }

  if (verbose) {
    message("Preprocessing complete")
  }

  return(object)
}

#' Find variable features
#'
#' @param object SpatialMetabolic object
#' @param n_features Number of features to select
#' @param method Method for finding variable features ("vst", "mvp", "disp")
#' @param verbose Print messages
#'
#' @return Object with variable features identified
#' @importFrom sparseMatrixStats rowVars rowMeans2
#' @importFrom SummarizedExperiment rowData rowData<-
#' @export
findVariableFeatures <- function(object,
                                 n_features = 2000,
                                 method = "vst",
                                 verbose = TRUE) {

  if (method == "vst") {
    # Variance stabilizing transformation approach
    counts_mat <- counts(object)
    log_counts <- logcounts(object)

    # Calculate mean and variance
    gene_means <- rowMeans2(counts_mat)
    gene_vars <- rowVars(log_counts)

    # Fit mean-variance relationship
    # Simple approach: use standardized variance
    not_const <- gene_vars > 0

    if (sum(not_const) < n_features) {
      n_features <- sum(not_const)
      if (verbose) {
        warning("Only ", n_features, " genes have non-zero variance")
      }
    }

    # Standardized variance (variance/mean)
    gene_disp <- rep(0, length(gene_means))
    gene_disp[not_const] <- gene_vars[not_const] / (gene_means[not_const] + 1)

    # Select top variable genes
    top_indices <- order(gene_disp, decreasing = TRUE)[1:n_features]

  } else if (method == "mvp") {
    # Mean-variance plot method
    log_counts <- logcounts(object)

    gene_means <- rowMeans2(log_counts)
    gene_vars <- rowVars(log_counts)

    # Fit loess curve
    fit <- loess(gene_vars ~ gene_means, span = 0.3)
    expected_vars <- predict(fit, gene_means)

    # Residual variance
    residual_vars <- gene_vars - expected_vars

    # Select top genes
    top_indices <- order(residual_vars, decreasing = TRUE)[1:n_features]

  } else if (method == "disp") {
    # Dispersion-based method
    counts_mat <- counts(object)

    gene_means <- rowMeans2(counts_mat)
    gene_vars <- rowVars(counts_mat)

    # Dispersion (overdispersion parameter)
    # Using method of moments estimator
    gene_disp <- (gene_vars - gene_means) / (gene_means^2 + 0.01)
    gene_disp[gene_disp < 0] <- 0

    # Select top genes
    top_indices <- order(gene_disp, decreasing = TRUE)[1:n_features]

  } else {
    stop("Method must be one of: vst, mvp, disp")
  }

  # Mark variable features
  variable_features <- rep(FALSE, nrow(object))
  variable_features[top_indices] <- TRUE

  # Add to rowData
  rd <- rowData(object)
  rd$variable_feature <- variable_features
  rd$mean_expression <- gene_means
  rd$variance <- gene_vars

  if (method == "vst" || method == "disp") {
    rd$dispersion <- gene_disp
  }

  rowData(object) <- rd

  if (verbose) {
    message("Selected ", n_features, " variable features")
  }

  return(object)
}

#' Run PCA on spatial data
#'
#' @param object SpatialMetabolic object
#' @param n_components Number of components
#' @param features Features to use (default: variable features)
#' @param center Whether to center data
#' @param scale Whether to scale data
#' @param verbose Print messages
#'
#' @return Object with PCA results
#' @importFrom stats prcomp
#' @importFrom SingleCellExperiment reducedDim<-
#' @export
runPCA <- function(object,
                   n_components = 30,
                   features = NULL,
                   center = TRUE,
                   scale = TRUE,
                   verbose = TRUE) {

  # Get features
  if (is.null(features)) {
    if ("variable_feature" %in% colnames(rowData(object))) {
      features <- rownames(object)[rowData(object)$variable_feature]
      if (verbose) {
        message("Using ", length(features), " variable features for PCA")
      }
    } else {
      stop("No variable features found. Run findVariableFeatures first.")
    }
  }

  # Get expression matrix
  expr_mat <- logcounts(object)[features, ]

  # Check dimensions
  max_components <- min(nrow(expr_mat), ncol(expr_mat)) - 1
  if (n_components > max_components) {
    n_components <- max_components
    if (verbose) {
      warning("Reducing n_components to ", n_components,
              " (max possible)")
    }
  }

  # Transpose for prcomp (expects samples in rows)
  expr_mat <- t(expr_mat)

  # Run PCA
  if (verbose) {
    message("Running PCA...")
  }

  pca_res <- prcomp(expr_mat,
                    center = center,
                    scale. = scale,
                    rank. = n_components)

  # Store results
  reducedDim(object, "PCA") <- pca_res$x

  # Store loadings and variance explained
  metadata(object)$PCA <- list(
    rotation = pca_res$rotation,
    sdev = pca_res$sdev,
    var_explained = pca_res$sdev^2 / sum(pca_res$sdev^2),
    features = features
  )

  if (verbose) {
    var_exp <- round(metadata(object)$PCA$var_explained[1:min(5, n_components)] * 100, 1)
    message("Top 5 PCs explain: ", paste(var_exp, "%", sep = "", collapse = ", "))
  }

  return(object)
}

#' Run UMAP on spatial data
#'
#' @param object SpatialMetabolic object
#' @param n_components Number of UMAP components
#' @param use_dims Number of PCs to use
#' @param n_neighbors Number of neighbors
#' @param min_dist Minimum distance
#' @param metric Distance metric
#' @param seed Random seed
#' @param verbose Print messages
#'
#' @return Object with UMAP results
#' @export
runUMAP <- function(object,
                    n_components = 2,
                    use_dims = 30,
                    n_neighbors = 30,
                    min_dist = 0.3,
                    metric = "euclidean",
                    seed = 42,
                    verbose = TRUE) {

  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' required for UMAP")
  }

  # Check for PCA
  if (!"PCA" %in% reducedDimNames(object)) {
    stop("PCA not found. Run runPCA first.")
  }

  # Get PCA embeddings
  pca_embed <- reducedDim(object, "PCA")
  use_dims <- min(use_dims, ncol(pca_embed))

  if (verbose) {
    message("Running UMAP on ", use_dims, " PCs...")
  }

  # Run UMAP
  set.seed(seed)
  umap_res <- uwot::umap(
    pca_embed[, 1:use_dims],
    n_components = n_components,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    verbose = verbose
  )

  # Store results
  colnames(umap_res) <- paste0("UMAP_", seq_len(n_components))
  reducedDim(object, "UMAP") <- umap_res

  return(object)
}

#' Subset SpatialMetabolic object
#'
#' @param object SpatialMetabolic object
#' @param cells Cells to keep
#' @param features Features to keep
#' @param samples Samples to keep
#'
#' @return Subsetted object
#' @export
subsetSpatial <- function(object,
                          cells = NULL,
                          features = NULL,
                          samples = NULL) {

  # Determine cells to keep
  if (!is.null(samples)) {
    sample_cells <- colData(object)$sample_id %in% samples
    if (!is.null(cells)) {
      cells <- intersect(cells, which(sample_cells))
    } else {
      cells <- which(sample_cells)
    }
  }

  # Subset
  if (!is.null(features) && !is.null(cells)) {
    object <- object[features, cells]
  } else if (!is.null(features)) {
    object <- object[features, ]
  } else if (!is.null(cells)) {
    object <- object[, cells]
  }

  # Update slots
  if (!is.null(cells)) {
    # Subset metabolic scores
    if (nrow(metabolicScores(object)) > 0) {
      metabolicScores(object) <- metabolicScores(object)[, cells, drop = FALSE]
    }

    # Subset comparison groups
    comparisonGroups(object) <- droplevels(comparisonGroups(object)[cells])
  }

  return(object)
}

#' Merge multiple SpatialMetabolic objects
#'
#' @param object_list List of SpatialMetabolic objects
#' @param merge_pathways How to handle pathways ("union", "intersection", "first")
#' @param recalculate_scores Whether to recalculate metabolic scores
#'
#' @return Merged SpatialMetabolic object
#' @export
mergeSpatial <- function(object_list,
                         merge_pathways = "union",
                         recalculate_scores = FALSE) {

  if (!is.list(object_list)) {
    stop("object_list must be a list of SpatialMetabolic objects")
  }

  if (length(object_list) < 2) {
    stop("Need at least 2 objects to merge")
  }

  # Check all are SpatialMetabolic
  if (!all(sapply(object_list, inherits, "SpatialMetabolic"))) {
    stop("All objects must be SpatialMetabolic objects")
  }

  # Use cbind from SpatialExperiment
  merged_spe <- do.call(cbind, object_list)

  # Convert back to SpatialMetabolic
  merged <- as(merged_spe, "SpatialMetabolic")

  # Handle pathways
  all_pathways <- lapply(object_list, metabolicPathways)

  if (merge_pathways == "union") {
    # Combine all unique pathways
    pathway_names <- unique(unlist(lapply(all_pathways, names)))
    merged_pathways <- list()

    for (pname in pathway_names) {
      # Get all gene sets for this pathway
      gene_sets <- lapply(all_pathways, function(x) x[[pname]])
      gene_sets <- gene_sets[!sapply(gene_sets, is.null)]

      # Union of genes
      merged_pathways[[pname]] <- unique(unlist(gene_sets))
    }

    metabolicPathways(merged) <- merged_pathways

  } else if (merge_pathways == "intersection") {
    # Only keep pathways present in all objects
    common_names <- Reduce(intersect, lapply(all_pathways, names))

    merged_pathways <- list()
    for (pname in common_names) {
      # Intersection of genes
      gene_sets <- lapply(all_pathways, function(x) x[[pname]])
      merged_pathways[[pname]] <- Reduce(intersect, gene_sets)
    }

    metabolicPathways(merged) <- merged_pathways

  } else if (merge_pathways == "first") {
    # Use pathways from first object
    metabolicPathways(merged) <- all_pathways[[1]]
  }

  # Recalculate scores if requested
  if (recalculate_scores && length(metabolicPathways(merged)) > 0) {
    merged <- calculateMetabolicScores(merged)
  }

  return(merged)
}

#' Export SpatialMetabolic to Seurat
#'
#' @param object SpatialMetabolic object
#' @param assay_name Name for the RNA assay
#' @param add_metabolic_scores Whether to add metabolic scores as metadata
#'
#' @return Seurat object
#' @export
toSeurat <- function(object,
                     assay_name = "RNA",
                     add_metabolic_scores = TRUE) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' required for conversion")
  }

  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts(object),
    meta.data = as.data.frame(colData(object)),
    assay = assay_name
  )

  # Add normalized data if present
  if (!all(logcounts(object) == 0)) {
    seurat_obj <- Seurat::SetAssayData(
      seurat_obj,
      slot = "data",
      new.data = logcounts(object)
    )
  }

  # Add spatial coordinates
  coords <- spatialCoords(object)
  colnames(coords) <- c("x", "y")

  # Create images list (simplified - would need proper image handling)
  image <- new("VisiumV1",
               coordinates = coords,
               spot.radius = 1)

  seurat_obj@images$slice1 <- image

  # Add metabolic scores as metadata
  if (add_metabolic_scores && nrow(metabolicScores(object)) > 0) {
    score_df <- as.data.frame(t(metabolicScores(object)))
    colnames(score_df) <- paste0("MetScore_", colnames(score_df))

    seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = score_df)
  }

  # Add dimension reductions
  if ("PCA" %in% reducedDimNames(object)) {
    seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
      embeddings = reducedDim(object, "PCA"),
      key = "PC_",
      assay = assay_name
    )
  }

  if ("UMAP" %in% reducedDimNames(object)) {
    seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = reducedDim(object, "UMAP"),
      key = "UMAP_",
      assay = assay_name
    )
  }

  return(seurat_obj)
}

#' Calculate quality control metrics
#'
#' @param object SpatialMetabolic object
#' @param mitochondrial_pattern Pattern to identify mitochondrial genes
#' @param ribosomal_pattern Pattern to identify ribosomal genes
#'
#' @return Object with QC metrics added to colData
#' @importFrom Matrix colSums
#' @export
calculateQCMetrics <- function(object,
                               mitochondrial_pattern = "^mt-|^Mt-|^MT-",
                               ribosomal_pattern = "^Rp[sl]|^RP[SL]") {

  counts_mat <- counts(object)

  # Basic metrics
  colData(object)$nCount_RNA <- colSums(counts_mat)
  colData(object)$nFeature_RNA <- colSums(counts_mat > 0)

  # Mitochondrial percentage
  mito_genes <- grep(mitochondrial_pattern, rownames(object), value = TRUE)
  if (length(mito_genes) > 0) {
    mito_counts <- colSums(counts_mat[mito_genes, , drop = FALSE])
    colData(object)$percent_mt <- mito_counts / colData(object)$nCount_RNA * 100
  } else {
    colData(object)$percent_mt <- 0
  }

  # Ribosomal percentage
  ribo_genes <- grep(ribosomal_pattern, rownames(object), value = TRUE)
  if (length(ribo_genes) > 0) {
    ribo_counts <- colSums(counts_mat[ribo_genes, , drop = FALSE])
    colData(object)$percent_ribo <- ribo_counts / colData(object)$nCount_RNA * 100
  } else {
    colData(object)$percent_ribo <- 0
  }

  # Log10 metrics for visualization
  colData(object)$log10_nCount_RNA <- log10(colData(object)$nCount_RNA + 1)
  colData(object)$log10_nFeature_RNA <- log10(colData(object)$nFeature_RNA + 1)

  message("Added QC metrics: nCount_RNA, nFeature_RNA, percent_mt, percent_ribo")

  return(object)
}

#' Save SpatialMetabolic object
#'
#' @param object SpatialMetabolic object
#' @param file Output file path
#' @param compress Whether to compress the file
#'
#' @export
saveSpatialMetabolic <- function(object, file, compress = TRUE) {
  if (compress) {
    saveRDS(object, file = file, compress = "gzip")
  } else {
    saveRDS(object, file = file, compress = FALSE)
  }

  message("Saved SpatialMetabolic object to: ", file)
}

#' Load SpatialMetabolic object
#'
#' @param file Input file path
#'
#' @return SpatialMetabolic object
#' @export
loadSpatialMetabolic <- function(file) {
  object <- readRDS(file)

  if (!inherits(object, "SpatialMetabolic")) {
    stop("File does not contain a SpatialMetabolic object")
  }

  return(object)
}
