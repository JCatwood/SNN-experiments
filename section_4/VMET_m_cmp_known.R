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
y_pred_VT <- y_obs
sd_pred_VT <- rep(0, length(y_obs))
time_VT <- 0
y_samp_VT <- matrix(y_obs,
                    nrow = length(y_obs),
                    ncol = n_samp, byrow = FALSE
)
for (m in m_seq) {
  cat("VMET sampling...\n")
  bgn_time <- Sys.time()
  for (i in 1:5) {
    for (j in 1:5) {
      cat("Partition", i, j, "\n")
      mask_ij_env <- locs[, 1] >= ((i - 1) * 0.2 - 0.015) &
        locs[, 1] <= (i * 0.2 + 0.015) &
        locs[, 2] >= ((j - 1) * 0.2 - 0.015) &
        locs[, 2] <= (j * 0.2 + 0.015)
      locs_ij_env <- locs[mask_ij_env, , drop = FALSE]
      y_obs_ij_env <- y_obs[mask_ij_env]
      mask_cens_ij_env <- mask_cens[mask_ij_env]
      cens_ub_ij_env <- cens_ub[mask_ij_env]
      cens_lb_ij_env <- cens_lb[mask_ij_env]
      covmat_ij_env <- covmat[mask_ij_env, mask_ij_env, drop = FALSE]
      locs_cens_ij_env <- locs_ij_env[mask_cens_ij_env, , drop = FALSE]
      mask_inner_ij_env <- locs_ij_env[, 1] >= (i - 1) * 0.2 &
        locs_ij_env[, 1] <= i * 0.2 &
        locs_ij_env[, 2] >= (j - 1) * 0.2 & locs_ij_env[, 2] <= j * 0.2
      set.seed(123)
      # if all responses are censored
      if (sum(mask_cens_ij_env) == length(y_obs_ij_env)) {
        samp_ij_env_VT <- VeccTMVN::mvrandn(
          lower = cens_lb_ij_env, upper = cens_ub_ij_env, mean = 0,
          sigma = covmat_ij_env,
          m = min(m, length(mask_cens_ij_env) - 1), N = n_samp
        )
      } else {
        tmp_mat <- solve(
          covmat_ij_env[!mask_cens_ij_env, !mask_cens_ij_env, drop = FALSE],
          covmat_ij_env[!mask_cens_ij_env, mask_cens_ij_env, drop = FALSE]
        )
        cond_covmat_ij_env_cens <-
          covmat_ij_env[mask_cens_ij_env, mask_cens_ij_env, drop = FALSE] -
          covmat_ij_env[mask_cens_ij_env, !mask_cens_ij_env, drop = FALSE] %*%
          tmp_mat
        cond_mean_ij_env_cens <- as.vector(
          t(y_obs_ij_env[!mask_cens_ij_env]) %*% tmp_mat
        )
        samp_ij_env_VT <- VeccTMVN::mvrandn(
          lower = cens_lb_ij_env[mask_cens_ij_env],
          upper = cens_ub_ij_env[mask_cens_ij_env],
          mean = cond_mean_ij_env_cens,
          sigma = cond_covmat_ij_env_cens,
          m = min(m, length(cond_mean_ij_env_cens) - 1), N = n_samp
        )
      }
      locs_cens_ij_env <- locs_ij_env[mask_cens_ij_env, , drop = FALSE]
      mask_cens_inner_ij_env <- locs_cens_ij_env[, 1] >= (i - 1) * 0.2 &
        locs_cens_ij_env[, 1] <= i * 0.2 &
        locs_cens_ij_env[, 2] >= (j - 1) * 0.2 & locs_cens_ij_env[, 2] <= j * 0.2
      y_pred_VT[mask_ij_env][mask_inner_ij_env & mask_cens_ij_env] <-
        rowMeans(samp_ij_env_VT)[mask_cens_inner_ij_env]
      sd_pred_VT[mask_ij_env][mask_inner_ij_env & mask_cens_ij_env] <-
        apply(samp_ij_env_VT, 1, sd)[mask_cens_inner_ij_env] / sqrt(n_samp)
      y_samp_VT[mask_ij_env, ][mask_inner_ij_env & mask_cens_ij_env, ] <-
        samp_ij_env_VT[mask_cens_inner_ij_env, ]
    }
  }
  end_time <- Sys.time()
  cat("VMET sampling is done\n")
  time_VT <- difftime(end_time, bgn_time, units = "secs")[[1]]
  cat(
    "> ", scene_ID, ",", m, ", RMSE, VT, known, ",
    sqrt(mean((y[mask_cens] - y_pred_VT[mask_cens])^2)), "\n"
  )
  cat(
    "> ", scene_ID, ",", m, ", NLL, VT, known, ",
    -mean(dnorm(y[mask_cens],
                mean = y_pred_VT[mask_cens],
                sd = sd_pred_VT[mask_cens]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ",", m, ", CRPS, VT, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens], dat = y_samp_VT[mask_cens, ]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ",", m, ", time, VT, known, ",
    time_VT, "\n"
  )
}

