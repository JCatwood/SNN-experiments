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
use_maxmin_order <- 0 # 0 no reorder, 1 maximin
use_snn_order <- 0
n_samp <- 50 # samples generated for posterior inference
run_seq_Vecc <- TRUE
use_parallel <- FALSE
args <- commandArgs(trailingOnly = TRUE)
# use command line args when running in batch on clusters
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2]) # simulation scenario ID
}

# data simulation ----------------------
source("data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- NA

# nntmvn ---------------------------------------
if (run_seq_Vecc) {
  if (use_maxmin_order == 0) {
    order <- 1:n
  } else if (use_maxmin_order == 1) {
    order <- GpGp::order_maxmin(locs)
  }
  y_obs_order <- y_obs[order]
  cens_ub_order <- cens_ub[order]
  cens_lb_order <- cens_lb[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  covmat_order <- covmat[order, order]
  for (m in m_seq) {
    cat("SNN sampling...\n")
    bgn_time <- Sys.time()
    if (use_parallel) {
      ncores <- 4
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      y_samp_seq_Vecc_order <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
        nntmvn::rptmvn(y_obs_order, cens_lb_order, cens_ub_order,
          mask_cens_order,
          m = m,
          covmat = covmat_order, locs = locs_order, ordering = use_snn_order,
          seed = i
        )
      }
      stopCluster(cl)
    } else {
      y_samp_seq_Vecc_order <- lapply(1:n_samp, function(seed_id) {
        nntmvn::rptmvn(y_obs_order, cens_lb_order, cens_ub_order,
          mask_cens_order,
          m = m,
          covmat = covmat_order, locs = locs_order, ordering = use_snn_order,
          seed = seed_id
        )
      })
    }
    end_time <- Sys.time()
    cat("SNN sampling is done\n")
    time_seq_Vecc <- difftime(end_time, bgn_time, units = "secs")[[1]]
    rev_order <- 1:n
    rev_order[order] <- 1:n
    y_samp_seq_Vecc <- matrix(unlist(y_samp_seq_Vecc_order), n,
      n_samp,
      byrow = FALSE
    )[rev_order, , drop = FALSE]
    y_pred_seq_Vecc <- rowMeans(y_samp_seq_Vecc)
    y_pred_cens_seq_Vecc <- y_pred_seq_Vecc[mask_cens]
    if (!file.exists("results")) {
      dir.create("results")
    }
    save(time_seq_Vecc, y_samp_seq_Vecc, file = paste0(
      "results/m_cmp_known_SNN_scene",
      scene_ID, "_m", m, "_order", use_snn_order, "_rep", k, ".RData"
    ))
    if (use_snn_order == 0) {
      cat(
        "> ", scene_ID, ",", m, ", RMSE, SNN, known, ",
        sqrt(mean((y[mask_cens] - y_pred_cens_seq_Vecc)^2)), "\n"
      )
      sd_cens_seq_Vecc <- apply(y_samp_seq_Vecc, 1, sd)[mask_cens] /
        sqrt(n_samp)
      cat(
        "> ", scene_ID, ",", m, ", NLL, SNN, known, ",
        -mean(dnorm(y[mask_cens],
          mean = y_pred_cens_seq_Vecc,
          sd = sd_cens_seq_Vecc
        )), "\n"
      )
      cat(
        "> ", scene_ID, ",", m, ", CRPS, SNN, known, ",
        mean(scoringRules::crps_sample(
          y = y[mask_cens],
          dat = y_samp_seq_Vecc[mask_cens, , drop = FALSE]
        )), "\n"
      )
    } else {
      cat(
        "> ", scene_ID, ",", m, ", RMSE, SNN_order, known, ",
        sqrt(mean((y[mask_cens] - y_pred_cens_seq_Vecc)^2)), "\n"
      )
      sd_cens_seq_Vecc <- apply(y_samp_seq_Vecc, 1, sd)[mask_cens] /
        sqrt(n_samp)
      cat(
        "> ", scene_ID, ",", m, ", NLL, SNN_order, known, ",
        -mean(dnorm(y[mask_cens],
          mean = y_pred_cens_seq_Vecc,
          sd = sd_cens_seq_Vecc
        )), "\n"
      )
      cat(
        "> ", scene_ID, ",", m, ", CRPS, SNN_order, known, ",
        mean(scoringRules::crps_sample(
          y = y[mask_cens],
          dat = y_samp_seq_Vecc[mask_cens, , drop = FALSE]
        )), "\n"
      )
    }
  }
}
