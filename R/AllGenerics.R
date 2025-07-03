#' @title Generic method definitions for SpatialMetabolic
#' @description This file contains all generic method definitions
#' @importFrom methods setGeneric

#' Access metabolic pathways
#' @param x SpatialMetabolic object
#' @export
setGeneric("metabolicPathways", function(x) standardGeneric("metabolicPathways"))

#' Set metabolic pathways
#' @param x SpatialMetabolic object
#' @param value List of pathways
#' @export
setGeneric("metabolicPathways<-", function(x, value) standardGeneric("metabolicPathways<-"))

#' Access metabolic scores
#' @param x SpatialMetabolic object
#' @export
setGeneric("metabolicScores", function(x) standardGeneric("metabolicScores"))

#' Set metabolic scores
#' @param x SpatialMetabolic object
#' @param value Matrix of scores
#' @export
setGeneric("metabolicScores<-", function(x, value) standardGeneric("metabolicScores<-"))

#' Access spatial graphs
#' @param x SpatialMetabolic object
#' @export
setGeneric("spatialGraphs", function(x) standardGeneric("spatialGraphs"))

#' Set spatial graphs
#' @param x SpatialMetabolic object
#' @param value List of graphs
#' @export
setGeneric("spatialGraphs<-", function(x, value) standardGeneric("spatialGraphs<-"))

#' Access analysis results
#' @param x SpatialMetabolic object
#' @export
setGeneric("analysisResults", function(x) standardGeneric("analysisResults"))

#' Set analysis results
#' @param x SpatialMetabolic object
#' @param value List of results
#' @export
setGeneric("analysisResults<-", function(x, value) standardGeneric("analysisResults<-"))

#' Access comparison groups
#' @param x SpatialMetabolic object
#' @export
setGeneric("comparisonGroups", function(x) standardGeneric("comparisonGroups"))

#' Set comparison groups
#' @param x SpatialMetabolic object
#' @param value Factor of groups
#' @export
setGeneric("comparisonGroups<-", function(x, value) standardGeneric("comparisonGroups<-"))

#' Calculate metabolic scores
#' @param object SpatialMetabolic object
#' @param ... Additional parameters
#' @export
setGeneric("calculateMetabolicScores", function(object, ...) standardGeneric("calculateMetabolicScores"))

#' Detect metabolic gradients
#' @param object SpatialMetabolic object
#' @param ... Additional parameters
#' @export
setGeneric("detectMetabolicGradients", function(object, ...) standardGeneric("detectMetabolicGradients"))

#' Compare metabolic states
#' @param object SpatialMetabolic object
#' @param ... Additional parameters
#' @export
setGeneric("compareMetabolicStates", function(object, ...) standardGeneric("compareMetabolicStates"))

#' Find metabolic neighborhoods
#' @param object SpatialMetabolic object
#' @param ... Additional parameters
#' @export
setGeneric("findMetabolicNeighborhoods", function(object, ...) standardGeneric("findMetabolicNeighborhoods"))
