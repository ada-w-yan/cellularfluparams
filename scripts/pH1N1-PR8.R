devtools::load_all("~/git_repos/cellularfluparams/")
library(lazymcmc)
library(parallel)
run_exp("pH1N1-PR8", "main", FALSE)