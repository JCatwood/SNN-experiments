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
use_snn_order <- 0
n_samp <- 50 # samples generated for posterior inference
run_SNN <- TRUE
use_parallel <- FALSE
args <- commandArgs(trailingOnly = TRUE)
# use command line args when running in batch on clusters
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2]) # simulation scenario ID
  use_snn_order <- as.integer(args[3]) # 0, 1, 2
}

# data simulation ----------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- NA
L <- t(chol(covmat))

source("../utils/score_output.R")

# nntmvn ---------------------------------------
if (run_SNN) {
  order <- 1:n
  y_obs_order <- y_obs[order]
  cens_ub_order <- cens_ub[order]
  cens_lb_order <- cens_lb[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  covmat_order <- covmat[order, order]
  for (m in m_seq) {
    cat("SNN sampling...\n")
    bgn_time <- Sys.time()
    NN <- if (scene_ID %in% c(2)) nntmvn::corr_nn(covmat_order, m - 1) else NULL
    if (use_parallel) {
      ncores <- 4
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      y_samp_SNN_order <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
        nntmvn::rptmvn(y_obs_order, cens_lb_order, cens_ub_order,
          mask_cens_order,
          m = m,
          covmat = covmat_order, locs = locs_order, NN = NN, ordering = use_snn_order,
          seed = i
        )
      }
      stopCluster(cl)
    } else {
      y_samp_SNN_order <- lapply(1:n_samp, function(seed_id) {
        nntmvn::rptmvn(y_obs_order, cens_lb_order, cens_ub_order,
          mask_cens_order,
          m = m,
          covmat = covmat_order, locs = locs_order, NN = NN, ordering = use_snn_order,
          seed = seed_id
        )
      })
    }
    end_time <- Sys.time()
    cat("SNN sampling is done\n")
    comp_time <- difftime(end_time, bgn_time, units = "secs")[[1]]
    rev_order <- 1:n
    rev_order[order] <- 1:n
    y_samp <- matrix(unlist(y_samp_SNN_order), n,
      n_samp,
      byrow = FALSE
    )[rev_order, , drop = FALSE]

    if (use_snn_order == 0) {
      method <- "SNN"
    } else if (use_snn_order == 1) {
      method <- "SNN_order_desc"
    } else {
      method <- "SNN_order_maximin"
    }
    
    score_output(
      y_samp[mask_cens, , drop = FALSE],
      y[mask_cens],
      comp_time,
      scene_ID = scene_ID,
      m = m,
      method = method,
      parms = "known"
    )
  }
}
