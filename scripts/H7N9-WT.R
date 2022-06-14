devtools::load_all("~/git_repos/cellularfluparams/")
library(lazymcmc)
library(parallel)
run_exp("H7N9-WT", "main", FALSE)