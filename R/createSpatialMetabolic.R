#' Create SpatialMetabolic object from Seurat objects
#'
#' @param seurat_list List of Seurat objects with spatial data
#' @param condition_col Column name containing condition information
#' @param sample_col Column name containing sample information
#' @param bin_size Bin size for Visium HD data (default: 8)
#' @param min_features Minimum number of features per spot
#' @param min_spots Minimum number of spots per feature
#'
#' @return SpatialMetabolic object
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment counts
#' @importFrom methods as
#' @export
#'
#' @examples
#' \dontrun{
#' # Load Seurat objects
#' wt1 <- Load10X_Spatial("data/WT1/")
#' wt2 <- Load10X_Spatial("data/WT2/")
#' svd1 <- Load10X_Spatial("data/SVD1/")
#' svd2 <- Load10X_Spatial("data/SVD2/")
#'
#' # Add metadata
#' wt1$condition <- "WT"
#' wt1$sample <- "WT1"
#' wt2$condition <- "WT"
#' wt2$sample <- "WT2"
#' svd1$condition <- "SVD"
#' svd1$sample <- "SVD1"
#' svd2$condition <- "SVD"
#' svd2$sample <- "SVD2"
#'
#' # Create SpatialMetabolic object
#' seurat_list <- list(wt1, wt2, svd1, svd2)
#' spm <- createSpatialMetabolic(seurat_list,
#'                               condition_col = "condition",
#'                               sample_col = "sample")
#' }
createSpatialMetabolic <- function(seurat_list,
                                   condition_col = "condition",
                                   sample_col = "sample",
                                   bin_size = 8,
                                   min_features = 200,
                                   min_spots = 10) {

  # Validate input
  if (!is.list(seurat_list)) {
    stop("seurat_list must be a list of Seurat objects")
  }

  if (length(seurat_list) == 0) {
    stop("seurat_list cannot be empty")
  }

  # Check that all are Seurat objects
  if (!all(sapply(seurat_list, function(x) inherits(x, "Seurat")))) {
    stop("All elements of seurat_list must be Seurat objects")
  }

  # Extract and combine data
  message("Extracting spatial data from ", length(seurat_list), " Seurat objects...")

  count_matrices <- list()
  spatial_coords <- list()
  cell_metadata <- list()

  for (i in seq_along(seurat_list)) {
    obj <- seurat_list[[i]]

    # Check for required metadata
    if (!condition_col %in% colnames(obj@meta.data)) {
      stop("Column '", condition_col, "' not found in Seurat object ", i)
    }

    if (!sample_col %in% colnames(obj@meta.data)) {
      stop("Column '", sample_col, "' not found in Seurat object ", i)
    }

    # Extract counts
    count_matrices[[i]] <- obj[["RNA"]]$counts

    # Extract spatial coordinates
    coords <- GetTissueCoordinates(obj)
    coords$cell_id <- paste0("sample", i, "_", rownames(coords))
    spatial_coords[[i]] <- coords

    # Extract metadata
    meta <- obj@meta.data
    meta$cell_id <- paste0("sample", i, "_", rownames(meta))
    meta$sample_id <- obj@meta.data[[sample_col]]
    meta$condition <- obj@meta.data[[condition_col]]
    cell_metadata[[i]] <- meta

    # Update column names for count matrix
    colnames(count_matrices[[i]]) <- paste0("sample", i, "_", colnames(count_matrices[[i]]))
  }

  message('combining all data')

  # Combine all data
  combined_counts <- do.call(cbind, count_matrices)
  combined_coords <- do.call(rbind, spatial_coords)
  combined_meta <- do.call(rbind, cell_metadata)

  message('verifying cell order')

  # Ensure matching order
  cell_order <- combined_coords$cell_id
  combined_meta <- combined_meta[match(cell_order, combined_meta$cell_id), ]
  combined_counts <- combined_counts[, match(cell_order, colnames(combined_counts))]

  message('building SpatialExperiment object')

  # Create SpatialExperiment base
  spe <- SpatialExperiment(
    assays = list(counts = combined_counts),
    colData = combined_meta,
    checkDimnames = FALSE,
    spatialCoords = as.matrix(combined_coords[, c("x", "y")]),
    sample_id = combined_meta$sample_id
  )

  message('converting to SpatialMetabolic')

  # Convert to SpatialMetabolic
  spm <- as(spe, "SpatialMetabolic")

  # Set comparison groups
  comparisonGroups(spm) <- factor(combined_meta$condition)

  # Quality control filtering
  spm <- .filterSpatialMetabolic(spm, min_features, min_spots)

  message("Created SpatialMetabolic object with ", ncol(spm), " spots from ",
          length(unique(colData(spm)$sample_id)), " samples")

  return(spm)
}

#' Helper function for filtering SpatialMetabolic objects
#'
#' @param spm SpatialMetabolic object
#' @param min_features Minimum features per spot
#' @param min_spots Minimum spots per feature
#' @return Filtered SpatialMetabolic object
#' @importFrom SingleCellExperiment counts
#' @keywords internal
.filterSpatialMetabolic <- function(spm, min_features, min_spots) {
  # Calculate QC metrics
  counts_mat <- counts(spm)

  # Filter spots with too few features
  n_features_per_spot <- colSums(counts_mat > 0)
  spot_qc <- n_features_per_spot >= min_features

  # Filter features present in too few spots
  n_spots_per_feature <- rowSums(counts_mat > 0)
  feature_qc <- n_spots_per_feature >= min_spots

  # Apply filtering
  spm <- spm[feature_qc, spot_qc]

  # Update comparison groups
  comparisonGroups(spm) <- droplevels(comparisonGroups(spm))

  message("Filtered to ", nrow(spm), " features and ", ncol(spm), " spots")
  message("  Removed ", sum(!feature_qc), " features with < ", min_spots, " spots")
  message("  Removed ", sum(!spot_qc), " spots with < ", min_features, " features")

  return(spm)
}

#' Convert SpatialExperiment to SpatialMetabolic
#'
#' @param spe SpatialExperiment object
#' @param condition_col Column with condition information
#' @return SpatialMetabolic object
#' @importFrom methods as
#' @importFrom SummarizedExperiment colData
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert existing SpatialExperiment
#' spe <- readRDS("spatial_experiment.rds")
#' spm <- spatialExperimentToMetabolic(spe, condition_col = "condition")
#' }
spatialExperimentToMetabolic <- function(spe, condition_col = "condition") {
  if (!inherits(spe, "SpatialExperiment")) {
    stop("Input must be a SpatialExperiment object")
  }

  # Convert
  spm <- as(spe, "SpatialMetabolic")

  # Set comparison groups if available
  if (condition_col %in% colnames(colData(spe))) {
    comparisonGroups(spm) <- factor(colData(spe)[[condition_col]])
  }

  return(spm)
}
