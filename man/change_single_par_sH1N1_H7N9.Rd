% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_single_par.R
\name{change_single_par_sH1N1_H7N9}
\alias{change_single_par_sH1N1_H7N9}
\title{calculate how the cellular infection parameters change when a model parameter
is changed from the value for sH1N1-WT to that for H7N9-WT}
\usage{
change_single_par_sH1N1_H7N9(rescale_par_name, sample_idxs)
}
\arguments{
\item{rescale_par_name}{character:
"baseline", "p_inf", "delta", "c_inf", "k1", or "beta_inf".  If "baseline",
don't change any parameters and just calculate the value of the cellular 
infection parmaeter for sH1N1.  Otherwise, change the value of the named
parameter to that for H7N9, and calculate the value of the cellular
infectin parameter.}

\item{sample_idxs}{a vector of iteration numbers of the MCMC chains.
Each row of sample_idxs corresponds to an sH1N1 parameter set and a 
H7N9 parameter set.  So for each row of sample_idxs, we get one value for
each cellular infection parameter.}
}
\value{
a named vector with elements "log10_R_0", "gen_time" and "r"
}
\description{
calculate how the cellular infection parameters change when a model parameter
is changed from the value for sH1N1-WT to that for H7N9-WT
}
