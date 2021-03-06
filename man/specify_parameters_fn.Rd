% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_for_inference.R
\name{specify_parameters_fn}
\alias{specify_parameters_fn}
\title{make parameter table}
\usage{
specify_parameters_fn(data_source, n_L = 10, n_I = 10,
  inc_data_logical = rep(T, 5), loss = TRUE,
  p_over_time = "constant")
}
\arguments{
\item{data_source}{character string: "WT", "PR8"}

\item{n_L}{numeric vector of length 1. number of latent stages}

\item{n_I}{numer vector of length 1.  number of infectious stages}

\item{inc_data_logical}{logical vector of length 5.  indicates which of
\code{c("sc_inf", "sc_tot", "mc_inf", "mc_tot", "mock")} we have data for}

\item{loss}{logical vector of length 1.  if TRUE, fit model with loss of virus
due to entry into target cells. if FALSE, fit model without loss of virus
due to entry into target cells.}

\item{p_over_time}{"constant" for constant rate of virus production over
infected cell lifetime, "linear" for linearly increasing rate of virus production
over infected cell lifetime}
}
\value{
data frame with columns
values: default parameter value
names: internal name for parameter
fixed:  if 0, fit parameter, otherwise fix parameter
lower_bound: lower bound of parameter for fitting
upper_bound: upper bound of parameter for fitting
steps: width of proposal distribution in parameter space
names_plot: character vector which which to label parameter distirbution plot
}
\description{
make parameter table
}
