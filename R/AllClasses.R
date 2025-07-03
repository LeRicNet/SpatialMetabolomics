#' @title SpatialMetabolic Class
#' @description S4 class for spatial metabolic analysis extending SpatialExperiment
#'
#' @slot metabolic_pathways list of metabolic pathway gene sets
#' @slot spatial_graphs list of spatial neighborhood graphs
#' @slot metabolic_scores matrix of pathway activity scores
#' @slot analysis_results list of analysis results
#' @slot comparison_groups factor defining comparison groups
#'
#' @export
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @importFrom methods setClass
setClass("SpatialMetabolic",
         contains = "SpatialExperiment",
         slots = c(
           metabolic_pathways = "list",
           spatial_graphs = "list",
           metabolic_scores = "matrix",
           analysis_results = "list",
           comparison_groups = "factor"
         ),
         prototype = list(
           metabolic_pathways = list(),
           spatial_graphs = list(),
           metabolic_scores = matrix(nrow = 0, ncol = 0),
           analysis_results = list(),
           comparison_groups = factor()
         )
)

#' Show method for SpatialMetabolic objects
#'
#' @param object A SpatialMetabolic object
#' @importFrom methods show
#' @export
setMethod("show", "SpatialMetabolic", function(object) {
  callNextMethod()
  cat("Metabolic pathways:", length(object@metabolic_pathways), "\n")
  if (length(object@metabolic_pathways) > 0) {
    cat("  ", paste(names(object@metabolic_pathways), collapse = ", "), "\n")
  }
  cat("Metabolic scores calculated:", nrow(object@metabolic_scores) > 0, "\n")
  if (nrow(object@metabolic_scores) > 0) {
    cat("  Dimensions:", nrow(object@metabolic_scores), "spots x",
        ncol(object@metabolic_scores), "pathways\n")
  }
  cat("Comparison groups:", nlevels(object@comparison_groups), "\n")
  if (nlevels(object@comparison_groups) > 0) {
    cat("  ", paste(levels(object@comparison_groups), collapse = ", "), "\n")
  }
  cat("Analysis results:", length(object@analysis_results), "\n")
  if (length(object@analysis_results) > 0) {
    cat("  ", paste(names(object@analysis_results), collapse = ", "), "\n")
  }
})
