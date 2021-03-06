% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_for_inference.R
\name{get_data_and_generate_posterior_fn}
\alias{get_data_and_generate_posterior_fn}
\title{read data for experiment, and define functions to calculate likelihood}
\usage{
get_data_and_generate_posterior_fn(parTab, strain = "sH1N1-WT",
  inc_data_logical, p_over_time = "constant")
}
\arguments{
\item{strain}{name of strain}

\item{inc_data_logical}{logical vector of length 5.  indicates which of
\code{c("sc_inf", "sc_tot", "mc_inf", "mc_tot", "mock")} we have data for}

\item{p_over_time}{"constant" for constant rate of virus production over
infected cell lifetime, "linear" for linearly increasing rate of virus production
over infected cell lifetime}
}
\value{
list with four elements:
data: list of data frames containing the data for each experiment
CREATE_POSTERIOR_FUNC: wrapper to create function to calculate
probability of data given parameters
CREATE_PRIOR_FUNC: wrapper to create function to calculate
probability of data given prior
}
\description{
read data for experiment, and define functions to calculate likelihood
}
