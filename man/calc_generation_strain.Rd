% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_generation_distribution.R
\name{calc_generation_strain}
\alias{calc_generation_strain}
\title{calculate the generation number distribution at the time of peak viral load
for the maximum likelihood parameters for a given strain}
\usage{
calc_generation_strain(strain, sensitivity)
}
\arguments{
\item{strain}{character: sH1N1, pH1N1, H5N1, H7N9 for the WT strains;
ENG195, Turkey05 for pH1N1-PR8 and H5N1-PR8 respectively}

\item{sensitivity}{character: choice of "main" (model in main paper),
"n_L_n_I" (changing the number of latent and infectious stages), 
"loss" (ignoring loss of virus due to entry into target cells),
"mc" (use multi-cycle data only),
"linear_p" (production rate of virions increases linearly with time)}
}
\value{
a data frame with three columns.  strain contains the name of the strain;
gen has generation numbers 1:20; and V has the proportion of infectious virions
belonging to that generation at the time of peak viral load.
}
\description{
calculate the generation number distribution at the time of peak viral load
for the maximum likelihood parameters for a given strain
}
