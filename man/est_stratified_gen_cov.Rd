% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genetic_covariance.R
\name{est_stratified_gen_cov}
\alias{est_stratified_gen_cov}
\title{Genetic covariance between traits using LDScore regression}
\usage{
est_stratified_gen_cov(
  paths,
  ld_score_annotation = "REF/1000G_Phase3_baselineLD_v2.2_ldscores/",
  sample_prev,
  population_prev,
  trait_names
)
}
\arguments{
\item{paths}{Paths to munged summary statistics}

\item{ld_score_annotation}{LDScores for annotation}

\item{sample_prev}{Sample prevalences for GWASs. NA for continuous}

\item{population_prev}{Population prevalences for traits. NA for continuous}

\item{trait_names}{Specify names for traits}
}
\value{
Stratified genetic covariance matrix estimated by GenomicSEM::s_ldsc
}
\description{
Genetic covariance between traits using LDScore regression
}
\examples{
NA
}
