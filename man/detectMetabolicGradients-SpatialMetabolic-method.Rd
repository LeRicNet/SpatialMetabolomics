% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialGradients.R
\name{detectMetabolicGradients,SpatialMetabolic-method}
\alias{detectMetabolicGradients,SpatialMetabolic-method}
\title{Detect spatial metabolic gradients}
\usage{
\S4method{detectMetabolicGradients}{SpatialMetabolic}(
  object,
  features = NULL,
  method = "moran",
  n_neighbors = 6,
  graph_type = "knn",
  radius = NULL,
  permutations = 99,
  BPPARAM = SerialParam()
)
}
\arguments{
\item{object}{SpatialMetabolic object}

\item{features}{Features to test (default: all metabolic pathways)}

\item{method}{Method for spatial detection ("moran", "geary", "SPARK", "variogram")}

\item{n_neighbors}{Number of neighbors for spatial graph}

\item{graph_type}{Type of spatial graph ("knn", "delaunay", "radius")}

\item{radius}{Radius for radius-based graph (only used if graph_type = "radius")}

\item{permutations}{Number of permutations for significance testing}

\item{BPPARAM}{BiocParallel parameters}
}
\value{
SpatialMetabolic object with gradient results
}
\description{
Detect spatial metabolic gradients
}
\examples{
\dontrun{
# Detect gradients in metabolic scores
spm <- detectMetabolicGradients(spm)

# Test specific features
spm <- detectMetabolicGradients(spm,
                               features = c("Ndufa1", "Atp5a1"),
                               method = "moran")

# Use different spatial graph
spm <- detectMetabolicGradients(spm,
                               method = "moran",
                               graph_type = "delaunay")
}
}
