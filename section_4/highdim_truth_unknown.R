library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(scoringRules)

# simulation settings ---------------
rm(list = ls())
set.seed(123)
scene_ID <- 1
k <- 1
m <- 30 # number of nearest neighbors
n_samp <- 50 # samples generated for posterior inference
run_CB <- TRUE
args <- commandArgs(trailingOnly = TRUE)
# use command line args when running in batch on clusters
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2]) # simulation scenario ID
}
# CensSpBayes
n_burn <- 20000
n_iter_MC <- 25000
thin <- 5

# data simulation ----------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- cens_ub[mask_cens] # CensSpBayes does not allow NA
if (!exists("cov_name")) {
  cov_name <- "matern15_isotropic"
}

source("../utils/score_output.R")

# CenSpBayes ------------------------------
if (run_CB) {
  library(CensSpBayes)
  bgn_time <- Sys.time()
  inla.mats <- create_inla_mats(
    S = locs, S.pred = locs[mask_cens, ],
    offset = c(0.01, 0.2),
    cutoff = 0.05,
    max.edge = c(0.01, 0.1)
  )
  X.obs <- matrix(1, nrow(locs), 1)
  X.pred <- matrix(1, sum(mask_cens), 1)
  cat("CB sampling begins...\n")
  y_samp_CB <- CensSpBayes::CensSpBayes(
    Y = y_obs, S = locs, X = X.obs,
    cutoff.Y = cens_ub,
    S.pred = locs[mask_cens, ], X.pred = X.pred,
    inla.mats = inla.mats,
    rho.init = 0.1, rho.upper = 5,
    iters = n_iter_MC, burn = n_burn, thin = thin, ret_samp = TRUE
  )
  cat("CB sampling done\n")
  end_time <- Sys.time()
  time_CB <- difftime(end_time, bgn_time, units = "secs")[[1]]

  score_output(y_samp_CB$Y.pred.samp, y[mask_cens], time_CB,
    scene_ID = scene_ID, method = "CB", parms = "unknown"
  )
}
