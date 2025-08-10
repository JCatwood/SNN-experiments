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
y_test <- y_test_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- cens_ub[mask_cens] # CensSpBayes does not allow NA
if (!exists("cov_name")) {
  cov_name <- "matern15_isotropic"
}
L <- t(chol(covmat))

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
  ret_obj <- CensSpBayes::CensSpBayes(
    Y = y_obs, S = locs, X = X.obs,
    cutoff.Y = cens_ub,
    S.pred = locs[mask_cens, ], X.pred = X.pred,
    inla.mats = inla.mats,
    rho.init = 0.1, rho.upper = 5,
    iters = n_iter_MC, burn = n_burn, thin = thin, ret_samp = TRUE
  )
  y_samp_CB <- matrix(y_obs,
    nrow = length(y_obs),
    ncol = ncol(ret_obj$Y.pred.samp), byrow = FALSE
  )
  y_samp_CB[mask_cens, ] <- ret_obj$Y.pred.samp
  cat("CB sampling done\n")
  end_time <- Sys.time()
  time_CB <- difftime(end_time, bgn_time, units = "secs")[[1]]

  kriging_score_output(y_samp_CB, y_test, time_CB,
    scene_ID = scene_ID, method = "CB", parms = "unknown"
  )
}

# save data for heatmap -------------------------
if (k == 1) {
  if (!file.exists("plots")) {
    dir.create("plots")
  }
  if (run_CB) {
    y_all <- rep(NA, length(ind_train) + length(ind_test))
    y_all[ind_train] <- y_samp_CB[, 1]
    y_all[ind_test] <- y_test
    write.table(y_all,
      file = paste0(
        "plots/samp_cmp_known_CB_scene",
        scene_ID, ".csv"
      ),
      row.names = F, col.names = F
    )
  }
}
