library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(scoringRules)

# simulation settings ------------------------------
rm(list = ls())
set.seed(123)
scene_ID <- 1
m <- 30 # number of nearest neighbors
reorder <- 0 # 0 no reorder, 1 maximin
n_samp <- 50 # samples generated for posterior inference
run_seq_Vecc <- TRUE
run_VT <- TRUE
run_TN <- TRUE
use_parallel <- FALSE
plot_heatmap <- FALSE
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
} else {
  k <- 1
}

# data simulation ----------------------
source("data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- NA

# nntmvn ---------------------------------------
if (run_seq_Vecc) {
  if (reorder == 0) {
    order <- 1:n
  } else if (reorder == 1) {
    order <- GpGp::order_maxmin(locs)
  }
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
    y_samp_seq_Vecc_order <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
      nntmvn::rtmvn_snn(y_obs_order, cens_lb_order, cens_ub_order,
        mask_cens_order,
        m = m,
        covmat = covmat_order, locs = locs_order,
        seed = i
      )
    }
    stopCluster(cl)
  } else {
    y_samp_seq_Vecc_order <- lapply(1:n_samp, function(seed_id) {
      nntmvn::rtmvn_snn(y_obs_order, cens_lb_order, cens_ub_order,
        mask_cens_order,
        m = m,
        covmat = covmat_order, locs = locs_order,
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
    "results/sim_data_seq_Vecc_scene",
    scene_ID, "_m", m, "_order", reorder, "_rep", k, ".RData"
  ))
  cat(
    "> ", scene_ID, ", RMSE, SNN, known, ",
    sqrt(mean((y[mask_cens] - y_pred_cens_seq_Vecc)^2)), "\n"
  )
  sd_cens_seq_Vecc <- apply(y_samp_seq_Vecc, 1, sd)[mask_cens] /
    sqrt(n_samp)
  cat(
    "> ", scene_ID, ", NLL, SNN, known, ",
    -mean(dnorm(y[mask_cens],
      mean = y_pred_cens_seq_Vecc,
      sd = sd_cens_seq_Vecc
    )), "\n"
  )
  cat(
    "> ", scene_ID, ", CRPS, SNN, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens],
      dat = y_samp_seq_Vecc[mask_cens, , drop = FALSE]
    )), "\n"
  )
}

# VeccTMVN and TruncatedNormal -----------------------
if (run_VT || run_TN) {
  if (run_VT) {
    y_pred_VT <- y_obs
    sd_pred_VT <- rep(0, length(y_obs))
    y_samp_VT <- matrix(y_obs,
      nrow = length(y_obs),
      ncol = n_samp, byrow = FALSE
    )
    time_VT <- 0
  }
  if (run_TN) {
    y_pred_TN <- y_obs
    sd_pred_TN <- rep(0, length(y_obs))
    y_samp_TN <- matrix(y_obs,
      nrow = length(y_obs),
      ncol = n_samp, byrow = FALSE
    )
    time_TN <- 0
  }
  for (i in 1:5) {
    for (j in 1:5) {
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
      if (run_VT) {
        bgn_time <- Sys.time()
        set.seed(123)
        cat("VT sampling at i =", i, "j =", j, "\n")
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
        end_time <- Sys.time()
        time_VT <- time_VT + difftime(end_time, bgn_time, units = "secs")[[1]]
      }

      if (run_TN) {
        bgn_time <- Sys.time()
        set.seed(123)
        cat("TN sampling at i =", i, "j =", j, "\n")
        covmat_ij_env <- covmat[mask_ij_env, mask_ij_env, drop = FALSE]
        if (sum(mask_cens_ij_env) == sum(mask_ij_env)) {
          cond_mean_ij_env_cens <- rep(0, sum(mask_ij_env))
          cond_covmat_ij_env_cens <- covmat_ij_env
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
        }
        # guarantee symmetry
        cond_covmat_ij_env_cens[lower.tri(cond_covmat_ij_env_cens)] <-
          t(cond_covmat_ij_env_cens)[lower.tri(cond_covmat_ij_env_cens)]
        samp_ij_env_TN <- t(TruncatedNormal::rtmvnorm(
          n_samp, cond_mean_ij_env_cens,
          cond_covmat_ij_env_cens, cens_lb_ij_env[mask_cens_ij_env],
          cens_ub_ij_env[mask_cens_ij_env]
        ))
        locs_cens_ij_env <- locs_ij_env[mask_cens_ij_env, , drop = FALSE]
        mask_cens_inner_ij_env <- locs_cens_ij_env[, 1] >= (i - 1) * 0.2 &
          locs_cens_ij_env[, 1] <= i * 0.2 &
          locs_cens_ij_env[, 2] >= (j - 1) * 0.2 & 
          locs_cens_ij_env[, 2] <= j * 0.2
        y_pred_TN[mask_ij_env][mask_inner_ij_env & mask_cens_ij_env] <-
          rowMeans(samp_ij_env_TN)[mask_cens_inner_ij_env]
        sd_pred_TN[mask_ij_env][mask_inner_ij_env & mask_cens_ij_env] <-
          apply(samp_ij_env_TN, 1, sd)[mask_cens_inner_ij_env] / sqrt(n_samp)
        y_samp_TN[mask_ij_env, ][mask_inner_ij_env & mask_cens_ij_env, ] <-
          samp_ij_env_TN[mask_cens_inner_ij_env, ]
        end_time <- Sys.time()
        time_TN <- time_TN + difftime(end_time, bgn_time, units = "secs")[[1]]
      }
    }
  }
  if (run_VT) {
    if (!file.exists("results")) {
      dir.create("results")
    }
    save(y_pred_VT, sd_pred_VT, y_samp_VT, time_VT,
      file = paste0("results/sim_data_VT_scene", scene_ID, "_rep", k, ".RData")
    )
    cat(
      "> ", scene_ID, ", RMSE, VT, known, ",
      sqrt(mean((y[mask_cens] - y_pred_VT[mask_cens])^2)), "\n"
    )
    cat(
      "> ", scene_ID, ", NLL, VT, known, ",
      -mean(dnorm(y[mask_cens],
        mean = y_pred_VT[mask_cens],
        sd = sd_pred_VT[mask_cens]
      )), "\n"
    )
    cat(
      "> ", scene_ID, ", CRPS, VT, known, ",
      mean(scoringRules::crps_sample(
        y = y[mask_cens], dat = y_samp_VT[mask_cens, ]
      )), "\n"
    )
  }
  if (run_TN) {
    if (!file.exists("results")) {
      dir.create("results")
    }
    save(y_pred_TN, sd_pred_TN, y_samp_TN, time_TN,
      file = paste0("results/sim_data_TN_scene", scene_ID, "_rep", k, ".RData")
    )
    cat(
      "> ", scene_ID, ", RMSE, TN, known, ",
      sqrt(mean((y[mask_cens] - y_pred_TN[mask_cens])^2)), "\n"
    )
    cat(
      "> ", scene_ID, ", NLL, TN, known, ",
      -mean(dnorm(y[mask_cens],
        mean = y_pred_TN[mask_cens],
        sd = sd_pred_TN[mask_cens]
      )), "\n"
    )
    cat(
      "> ", scene_ID, ", CRPS, TN, known, ",
      mean(scoringRules::crps_sample(
        y = y[mask_cens], dat = y_samp_TN[mask_cens, ]
      )), "\n"
    )
  }
}

# heatmap -------------------------
if (plot_heatmap) {
  library(fields)
  load(paste0(
    "results/sim_data_seq_Vecc_scene",
    scene_ID, "_m", m, "_order", reorder, "_rep", k, ".RData"
  ))
  load(paste0("results/sim_data_VT_scene", scene_ID, "_rep", k, ".RData"))
  load(paste0("results/sim_data_TN_scene", scene_ID, "_rep", k, ".RData"))
  zlim <- range(y_samp_seq_Vecc, y_samp_TN, y_samp_VT)
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
      "plots/sim_data_seq_Vecc_scene", scene_ID, "_m", m,
      "_order", reorder, ".pdf"
    ),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_seq_Vecc, sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0("plots/sim_data_VT_scene", scene_ID, ".pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_VT, sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0("plots/sim_data_TN_scene", scene_ID, ".pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_TN, sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
}
