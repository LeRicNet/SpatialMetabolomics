% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createSpatialMetabolic.R
\name{spatialExperimentToMetabolic}
\alias{spatialExperimentToMetabolic}
\title{Convert SpatialExperiment to SpatialMetabolic}
\usage{
spatialExperimentToMetabolic(spe, condition_col = "condition")
}
\arguments{
\item{spe}{SpatialExperiment object}

\item{condition_col}{Column with condition information}
}
\value{
SpatialMetabolic object
}
\description{
Convert SpatialExperiment to SpatialMetabolic
}
\examples{
\dontrun{
# Convert existing SpatialExperiment
spe <- readRDS("spatial_experiment.rds")
spm <- spatialExperimentToMetabolic(spe, condition_col = "condition")
}
}
