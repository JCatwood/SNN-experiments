library(GpGp)
library(VeccTMVN)

rm(list = ls())
scene_ID <- 1
m <- 30
n_samp <- 50
subset_size <- 2500
n_burn <- 20000
n_iter_MC <- 25000
thin <- 5

args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args)) as.integer(args[1]) else 1

if (!requireNamespace("CensSpBayes", quietly = TRUE)) {
  stop("This script requires the CensSpBayes package.")
}

# Common simulated data ---------------------------------------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs_SNN <- y
y_obs_SNN[mask_cens] <- NA
y_obs_CSB <- y
y_obs_CSB[mask_cens] <- cens_ub[mask_cens]
source("../utils/score_output.R")

# CSB ---------------------------------------------------------------------
set.seed(123)
bgn_time <- Sys.time()
inla_mats <- CensSpBayes::create_inla_mats(
  S = locs, S.pred = locs[mask_cens, , drop = FALSE],
  offset = c(0.01, 0.2), cutoff = 0.05, max.edge = c(0.01, 0.1)
)
ret_CSB <- CensSpBayes::CensSpBayes(
  Y = y_obs_CSB, S = locs, X = matrix(1, n, 1),
  cutoff.Y = cens_ub,
  S.pred = locs[mask_cens, , drop = FALSE],
  X.pred = matrix(1, sum(mask_cens), 1),
  inla.mats = inla_mats, rho.init = 0.1, rho.upper = 5,
  iters = n_iter_MC, burn = n_burn, thin = thin, ret_samp = TRUE
)
comp_time_CSB <- difftime(Sys.time(), bgn_time, units = "secs")[[1]]
y_cens_samp_CSB <- tail(ret_CSB$Y.pred.samp, sum(mask_cens))
score_output(
  y_cens_samp_CSB,
  y[mask_cens],
  comp_time_CSB,
  scene_ID = scene_ID,
  method = "CSB",
  parms = "unknown"
)

