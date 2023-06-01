#######################################################################
##* Run BayesLSA on Application case of United States maps
##* 1. read the raw data
##* 2. run BayesLASA
#######################################################################

library(dplyr)
library(mcclust)
library(Rcpp)
library(doParallel)
library(foreach)
registerDoParallel(4)

#***************
#* file path
#***************
input <- "manuscript_reproducibility/results/real_data_us_state/"
input_raw <- "data/real_data_us_state/"
code_path <- "code/"

#* raw data
load(file.path(input_raw, "states_outline.Rdata"), verbose = T)
state_list <- names(outline_polygon)
state_pns <- unlist(lapply(outline_polygon, nrow))

# Identify landmarks using BayesLASA --------------------------------------

####################################################################################
##** NOTE: It takes a while to run the MCMC and identify the credible interval
##** Ready-to-use results saved in the 'manuscript_reproducibility/results/real_data_us_state/BayesLASA'
##*
####################################################################################
source(file.path(code_path, "landmark_detection/sim_polygon_gaussian.R"))
sourceCpp(file.path(code_path, "landmark_detection/MCMC_shape.cpp"))

load(file.path(input_raw, "states_outline.Rdata"), verbose = T)
state_list <- names(outline_polygon)
state_pns <- unlist(lapply(outline_polygon, nrow))
for (mm in 1:length(state_list)) {
  state_name <- state_list[mm]
  original_dat <- outline_polygon[[mm]]
  normalized <- pc_normalizor(original_dat)

  dat <- normalized$pc
  dat <- dat[-nrow(dat), ]
  n <- nrow(dat)
  beta_sigma <- c(0.01, 0.001, 0.0001, 0.00001)
  names(beta_sigma) <- c("0.01", "0.001", "0.0001", "0.00001")
  est.K <- 4
  fold <- max(5 * 1e4, 100 * n)

  hpts <- chull(dat)
  hpts <- sort(hpts)
  gamma_i <- rep(0, n)
  gamma_i[hpts] <- 1
  # gamma_i = L2gamma(hpts,n)

  foreach(k = 1:4) %dopar% {
    r <- MCMC_shape(dat,
      iter = fold, estK = est.K, gamma_i = gamma_i,
      alpha_sigma = 3, beta_sigma = beta_sigma[k], ppm_store = T
    )
    saveRDS(r, file = file.path(input, "BayesLASA", paste0("state_", mm, "_bs_", names(beta_sigma)[k], "_MCMCres.rds")))
    rm(r)
  }
}
