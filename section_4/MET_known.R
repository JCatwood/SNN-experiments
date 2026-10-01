library(TruncatedNormal)

# simulation settings ------------------------------
rm(list = ls())
n_samp <- 50
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  k <- as.integer(args[1])
  scene_ID <- as.integer(args[2])
} else {
  k <- 1
  scene_ID <- 1
}

# data simulation ----------------------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- NA

source("../utils/score_output.R")
source("sample_func_VT_TN.R")

# MET ----------------------------------------------
y_samp_MET <- matrix(y_obs, length(y_obs), n_samp)
bgn_time <- Sys.time()
set.seed(123)
ret_obj <- sample_wrapper(0, 1, 0, 1, "TN")
y_samp_MET[ret_obj$ind, ] <- ret_obj$samp
comp_time <- difftime(Sys.time(), bgn_time, units = "secs")[[1]]

score_output(
  y_samp_MET[mask_cens, , drop = FALSE],
  y[mask_cens],
  comp_time,
  scene_ID = scene_ID,
  method = "MET",
  parms = "known"
)
