% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_run_all.R
\name{run_exp}
\alias{run_exp}
\title{estimate parameters for the different models in the study}
\usage{
run_exp(strain, sensitivity, short)
}
\arguments{
\item{strain}{character: sH1N1, pH1N1, H5N1, H7N9 for the WT strains;
ENG195, Turkey05 for pH1N1-PR8 and H5N1-PR8 respectively}

\item{sensitivity}{character: choice of "main" (model in main paper),
"n_L_n_I" (changing the number of latent and infectious stages), 
"loss" (ignoring loss of virus due to entry into target cells),
"mc" (use multi-cycle data only),
"linear_p" (production rate of virions increases linearly with time)}

\item{short}{logical.  if TRUE, run short MCMC chain for testing purposes;
otherwise, run until convergence}
}
\value{
NULL (results saved to file)
}
\description{
estimate parameters for the different models in the study
}
