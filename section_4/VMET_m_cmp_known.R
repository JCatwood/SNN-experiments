library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(scoringRules)

# simulation settings ------------------------------
rm(list = ls())
set.seed(123)
scene_ID <- 1
k <- 1
m_seq <- seq(from = 10, to = 50, by = 10) # number of nearest neighbors
n_samp <- 50 # samples generated for posterior inference
args <- commandArgs(trailingOnly = TRUE)
# use command line args when running in batch on clusters
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2]) # simulation scenario ID
}

# data simulation ----------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- NA

# VMET ---------------------------------------
source("sample_func_VT_TN.R")
for (m in m_seq) {
  y_samp_VT <- matrix(y_obs,
    nrow = length(y_obs),
    ncol = n_samp, byrow = FALSE
  )
  bgn_time <- Sys.time()
  set.seed(123)
  ret_obj <- sample_wrapper(0, 1, 0, 1, "VT")
  y_samp_VT[ret_obj$ind, ] <- ret_obj$samp
  end_time <- Sys.time()
  time_VT <- difftime(end_time, bgn_time, units = "secs")[[1]]
  y_pred_VT <- rowMeans(y_samp_VT)
  cat(
    "> ", scene_ID, ",", m, ", RMSE, VT, known, ",
    sqrt(mean((y[mask_cens] - y_pred_VT[mask_cens])^2)), "\n"
  )
  cat(
    "> ", scene_ID, ",", m, ", CRPS, VT, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens], dat = y_samp_VT[mask_cens, , drop = FALSE]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ",", m, ", time, VT, known, ",
    time_VT, "\n"
  )
}
