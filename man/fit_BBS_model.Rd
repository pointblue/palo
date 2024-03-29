% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_BBS_model.R
\name{fit_BBS_model}
\alias{fit_BBS_model}
\title{Fit BBS-style hierarchical model for estimating trends in point count survey data}
\usage{
fit_BBS_model(inputdata, n.adapt = 1000, n.burnin = 5000,
  n.sample = 10000, n.chains = 3, overdispersion = TRUE, ...)
}
\arguments{
\item{inputdata}{Inputdata created by running \code{\link{setup_BBS_model}}}

\item{n.adapt}{Number of iterations for adaptation, defaults to 1000}

\item{n.burnin}{Number of iterations for burn-in, defaults to 5000}

\item{n.sample}{Number of iterations for sampling, defauts to 10000}

\item{n.chains}{Number of MCMC chains, defaults to 3}

\item{overdispersion}{Defaults to TRUE, to include a term in the model for
fitting overdispersion, as in BBS models}

\item{...}{Additional arguments passed to \code{\link[rjags]{jags.model}},
e.g. initial values.}
}
\value{
Returns a coda.samples object from rjags. Also displays summary
statistics from MCMCsummary and prints trace plots from MCMCplot.
}
\description{
Fit BBS-style hierarchical model for estimating trends in point count survey data
}
