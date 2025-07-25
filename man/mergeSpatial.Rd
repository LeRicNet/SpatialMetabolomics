% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mergeSpatial}
\alias{mergeSpatial}
\title{Merge multiple SpatialMetabolic objects}
\usage{
mergeSpatial(object_list, merge_pathways = "union", recalculate_scores = FALSE)
}
\arguments{
\item{object_list}{List of SpatialMetabolic objects}

\item{merge_pathways}{How to handle pathways ("union", "intersection", "first")}

\item{recalculate_scores}{Whether to recalculate metabolic scores}
}
\value{
Merged SpatialMetabolic object
}
\description{
Merge multiple SpatialMetabolic objects
}
