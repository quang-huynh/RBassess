% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{simulation}
\alias{simulation}
\title{Generates stochastic length and age-length samples to run a simulation}
\usage{
simulation(RBdata, pars, N = c(50, 100), CV_L_stock = 0,
  CV_stock_d = 0, seed = 1L)
}
\arguments{
\item{RBdata}{\linkS4class{RBdata} object used for setting up the simulation.}

\item{pars}{Matrix of parameters. See details.}

\item{N}{Numeric vector of length two for the number of age-length and length samples.}

\item{CV_L_stock}{Numeric, the CV of length at stocking among simulations.}

\item{CV_stock_d}{Numeric, the CV of stocking density among simulations.}

\item{seed}{Integer, for replicating random number generation.}
}
\description{
Simulate data for the assessment model.
}
\details{
The \code{pars} matrix should be nsim rows long and have columns.
}
