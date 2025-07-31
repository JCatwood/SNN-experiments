library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(scoringRules)

# simulation settings ------------------------------
rm(list = ls())
set.seed(123)
m <- 30 # number of nearest neighbors
n_samp <- 50 # samples generated for posterior inference
run_SNN <- TRUE
run_VT <- TRUE
run_TN <- TRUE
use_parallel <- FALSE
plot_heatmap <- FALSE
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2])
  use_snn_order <- as.integer(args[3]) # 0, 1, 2
} else {
  k <- 1
  scene_ID <- 1
  use_snn_order <- 0
}

# data simulation ----------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- NA

# nntmvn ---------------------------------------
if (run_SNN) {
  order <- 1:n
  y_obs_order <- y_obs[order]
  cens_ub_order <- cens_ub[order]
  cens_lb_order <- cens_lb[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  covmat_order <- covmat[order, order]
  cat("SNN sampling...\n")
  bgn_time <- Sys.time()
  if (use_parallel) {
    ncores <- 4
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    y_samp_SNN_order <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
      nntmvn::rptmvn(y_obs_order, cens_lb_order, cens_ub_order,
        mask_cens_order,
        m = m,
        covmat = covmat_order, locs = locs_order, ordering = use_snn_order,
        seed = i
      )
    }
    stopCluster(cl)
  } else {
    y_samp_SNN_order <- lapply(1:n_samp, function(seed_id) {
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
  time_SNN <- difftime(end_time, bgn_time, units = "secs")[[1]]
  rev_order <- 1:n
  rev_order[order] <- 1:n
  y_samp_SNN <- matrix(unlist(y_samp_SNN_order), n,
    n_samp,
    byrow = FALSE
  )[rev_order, , drop = FALSE]
  y_pred_SNN <- rowMeans(y_samp_SNN)
  y_pred_cens_SNN <- y_pred_SNN[mask_cens]
  cat(
    "> ", scene_ID, ", RMSE, SNN, known, ",
    sqrt(mean((y[mask_cens] - y_pred_cens_SNN)^2)), "\n"
  )
  cat(
    "> ", scene_ID, ", CRPS, SNN, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens],
      dat = y_samp_SNN[mask_cens, , drop = FALSE]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ", time, SNN, known, ",
    time_SNN, "\n"
  )
}

# VeccTMVN and TruncatedNormal -----------------------
if (run_VT || run_TN) {
  source("sample_func_VT_TN.R")
}

if (run_VT) {
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
    "> ", scene_ID, ", RMSE, VT, known, ",
    sqrt(mean((y[mask_cens] - y_pred_VT[mask_cens])^2)), "\n"
  )
  cat(
    "> ", scene_ID, ", CRPS, VT, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens], dat = y_samp_VT[mask_cens, , drop = FALSE]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ", time, VT, known, ",
    time_VT, "\n"
  )
}

if (run_TN) {
  y_samp_TN <- matrix(y_obs,
    nrow = length(y_obs),
    ncol = n_samp, byrow = FALSE
  )
  bgn_time <- Sys.time()
  set.seed(123)
  ret_obj <- sample_wrapper(0, 1, 0, 1, "TN")
  y_samp_TN[ret_obj$ind, ] <- ret_obj$samp
  end_time <- Sys.time()
  time_TN <- difftime(end_time, bgn_time, units = "secs")[[1]]
  y_pred_TN <- rowMeans(y_samp_TN)
  cat(
    "> ", scene_ID, ", RMSE, TN, known, ",
    sqrt(mean((y[mask_cens] - y_pred_TN[mask_cens])^2)), "\n"
  )
  cat(
    "> ", scene_ID, ", CRPS, TN, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens], dat = y_samp_TN[mask_cens, , drop = FALSE]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ", time, TN, known, ",
    time_TN, "\n"
  )
}

# heatmap -------------------------
if (plot_heatmap) {
  library(fields)
  zlim <- range(y_samp_SNN, y_samp_TN, y_samp_VT)
  if (!file.exists("plots")) {
    dir.create("plots")
  }
  pdf(
    file = paste0("plots/sim_data_scene", scene_ID, "_true.pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y, sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0(
      "plots/samp_cmp_known_SNN_scene", scene_ID, "_m", m,
      "_order", use_snn_order, ".pdf"
    ),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_SNN[, 1], sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0("plots/samp_cmp_known_VT_scene", scene_ID, ".pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_VT[, 1], sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0("plots/samp_cmp_known_TN_scene", scene_ID, ".pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_TN[, 1], sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
}
