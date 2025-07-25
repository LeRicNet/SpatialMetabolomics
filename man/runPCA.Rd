% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{runPCA}
\alias{runPCA}
\title{Run PCA on spatial data}
\usage{
runPCA(
  object,
  n_components = 30,
  features = NULL,
  center = TRUE,
  scale = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{object}{SpatialMetabolic object}

\item{n_components}{Number of components}

\item{features}{Features to use (default: variable features)}

\item{center}{Whether to center data}

\item{scale}{Whether to scale data}

\item{verbose}{Print messages}
}
\value{
Object with PCA results
}
\description{
Run PCA on spatial data
}
