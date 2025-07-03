#' @title Accessor methods for SpatialMetabolic
#' @description Implementation of accessor methods
#' @importFrom methods setMethod

#' @describeIn metabolicPathways Access metabolic pathways
#' @param x A SpatialMetabolic object
#' @return List of pathway gene sets
#' @export
setMethod("metabolicPathways", "SpatialMetabolic", function(x) {
  x@metabolic_pathways
})

#' @describeIn metabolicPathways Set metabolic pathways
#' @param x A SpatialMetabolic object
#' @param value List of pathway gene sets
#' @return Modified SpatialMetabolic object
#' @export
setMethod("metabolicPathways<-", "SpatialMetabolic", function(x, value) {
  if (!is.list(value)) {
    stop("value must be a list")
  }
  x@metabolic_pathways <- value
  validObject(x)
  x
})

#' @describeIn metabolicScores Access metabolic scores
#' @param x A SpatialMetabolic object
#' @return Matrix of metabolic scores
#' @export
setMethod("metabolicScores", "SpatialMetabolic", function(x) {
  x@metabolic_scores
})

#' @describeIn metabolicScores Set metabolic scores
#' @param x A SpatialMetabolic object
#' @param value Matrix of metabolic scores
#' @return Modified SpatialMetabolic object
#' @export
setMethod("metabolicScores<-", "SpatialMetabolic", function(x, value) {
  if (!is.matrix(value)) {
    stop("value must be a matrix")
  }
  x@metabolic_scores <- value
  validObject(x)
  x
})

#' @describeIn spatialGraphs Access spatial graphs
#' @param x A SpatialMetabolic object
#' @return List of spatial graphs
#' @export
setMethod("spatialGraphs", "SpatialMetabolic", function(x) {
  x@spatial_graphs
})

#' @describeIn spatialGraphs Set spatial graphs
#' @param x A SpatialMetabolic object
#' @param value List of spatial graphs
#' @return Modified SpatialMetabolic object
#' @export
setMethod("spatialGraphs<-", "SpatialMetabolic", function(x, value) {
  if (!is.list(value)) {
    stop("value must be a list")
  }
  x@spatial_graphs <- value
  validObject(x)
  x
})

#' @describeIn analysisResults Access analysis results
#' @param x A SpatialMetabolic object
#' @return List of analysis results
#' @export
setMethod("analysisResults", "SpatialMetabolic", function(x) {
  x@analysis_results
})

#' @describeIn analysisResults Set analysis results
#' @param x A SpatialMetabolic object
#' @param value List of analysis results
#' @return Modified SpatialMetabolic object
#' @export
setMethod("analysisResults<-", "SpatialMetabolic", function(x, value) {
  if (!is.list(value)) {
    stop("value must be a list")
  }
  x@analysis_results <- value
  validObject(x)
  x
})

#' @describeIn comparisonGroups Access comparison groups
#' @param x A SpatialMetabolic object
#' @return Factor of comparison groups
#' @export
setMethod("comparisonGroups", "SpatialMetabolic", function(x) {
  x@comparison_groups
})

#' @describeIn comparisonGroups Set comparison groups
#' @param x A SpatialMetabolic object
#' @param value Factor of comparison groups
#' @return Modified SpatialMetabolic object
#' @export
setMethod("comparisonGroups<-", "SpatialMetabolic", function(x, value) {
  if (!is.factor(value)) {
    value <- factor(value)
  }
  x@comparison_groups <- value
  validObject(x)
  x
})
