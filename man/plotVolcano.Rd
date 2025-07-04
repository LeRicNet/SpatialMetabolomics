% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotVolcano}
\alias{plotVolcano}
\title{Volcano plot for differential expression results}
\usage{
plotVolcano(
  de_results,
  pval_cutoff = 0.05,
  fc_cutoff = 0.5,
  label_top = 10,
  colors = c("blue", "grey", "red")
)
}
\arguments{
\item{de_results}{Results from compareMetabolicStates}

\item{pval_cutoff}{P-value cutoff}

\item{fc_cutoff}{Fold change cutoff (absolute value)}

\item{label_top}{Number of top features to label}

\item{colors}{Colors for down, non-sig, and up}
}
\value{
ggplot object
}
\description{
Volcano plot for differential expression results
}
