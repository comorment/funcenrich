% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genetic_correlation.R
\name{est_genetic_correlation}
\alias{est_genetic_correlation}
\title{Genetic correlation}
\usage{
est_genetic_correlation(trait1, trait2, gen_cov, adjust = NULL)
}
\arguments{
\item{trait1}{Name of trait 1}

\item{trait2}{Name of trait 2}

\item{gen_cov}{Genetic covariance matrix estimated by ldsc}

\item{adjust}{Name of trait to adjust for}
}
\value{
Tibble with genetic correlation between trait1 and trait2 adjusted for genetic components of trait in adjust
}
\description{
Genetic correlation
}
\examples{
NA
}
