library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(scoringRules)

# simulation settings ------------------------------
rm(list = ls())
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
# general
scene_ID <- 1 # Matern 1.5 kernel
m <- 30 # number of nearest neighbors
reorder <- 0 # 0 no reorder, 1 maximin
n_samp <- 50 # samples generated for posterior inference
run_seq_Vecc <- TRUE
run_VT <- TRUE
run_TN <- TRUE
use_parallel <- FALSE
plot_heatmap <- FALSE


# data simulation ----------------------
if (scene_ID == 1) {
  tmp_vec <- seq(from = 0, to = 1, length.out = 100)
  locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
  cov_func <- GpGp::matern15_isotropic
  cov_parms <- c(1.0, 0.03, 0.0001)
  cov_name <- "matern15_isotropic"
  covmat <- cov_func(cov_parms, locs)
  cat("Generating GP ...", "\n")
  y <- as.vector(t(chol(covmat)) %*% rnorm(nrow(locs)))
  cat("GP generated", "\n")
  n <- nrow(locs)
  levl_cens <- rep(1, n)
  rm(tmp_vec)
}

mask_cens <- y < levl_cens
y_obs <- y
y_obs[mask_cens] <- levl_cens[mask_cens]

# nntmvn ---------------------------------------
if (run_seq_Vecc) {
  if (reorder == 0) {
    order <- 1:n
  } else if (reorder == 1) {
    order <- GpGp::order_maxmin(locs)
  }
  y_obs_order <- y_obs[order]
  levl_cens_order <- levl_cens[order]
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
      nntmvn::rtmvn_snn(y_obs_order, rep(-Inf, n), levl_cens_order,
        mask_cens_order,
        m = m,
        covmat = covmat_order, locs = locs_order,
        seed = i
      )
    }
    stopCluster(cl)
  } else {
    y_samp_seq_Vecc_order <- lapply(1:n_samp, function(seed_id) {
      nntmvn::rtmvn_snn(y_obs_order, rep(-Inf, n), levl_cens_order,
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
    scene_ID, "_m", m, "_order", reorder, ".RData"
  ))
  cat(
    "RMSE for seq Vecc is ",
    sqrt(mean((y[mask_cens] - y_pred_cens_seq_Vecc)^2)), "\n"
  )
  sd_cens_seq_Vecc <- apply(y_samp_seq_Vecc, 1, sd)[mask_cens] /
    sqrt(n_samp)
  cat(
    "NLL for seq Vecc is ",
    -mean(dnorm(y[mask_cens],
      mean = y_pred_cens_seq_Vecc,
      sd = sd_cens_seq_Vecc
    )), "\n"
  )
  cat(
    "CRPS for seq Vecc is ",
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
      levl_cens_ij_env <- levl_cens[mask_ij_env]
      mask_inner_ij_env <- locs_ij_env[, 1] >= (i - 1) * 0.2 &
        locs_ij_env[, 1] <= i * 0.2 &
        locs_ij_env[, 2] >= (j - 1) * 0.2 & locs_ij_env[, 2] <= j * 0.2
      if (run_VT) {
        bgn_time <- Sys.time()
        set.seed(123)
        cat("VT sampling at i =", i, "j =", j, "\n")
        if (sum(mask_cens_ij_env) == length(y_obs_ij_env)) {
          samp_ij_env_VT <- VeccTMVN::mvrandn(
            lower = -Inf, upper = levl_cens_ij_env, mean = 0, locs = locs_ij_env,
            covName = cov_name, covParms = cov_parms,
            m = min(m, length(mask_cens_ij_env) - 1), N = n_samp
          )
        } else {
          samp_ij_env_VT <- VeccTMVN::ptmvrandn(
            locs_ij_env, which(mask_cens_ij_env),
            y_obs_ij_env, levl_cens_ij_env,
            cov_name, cov_parms,
            m = min(m, length(mask_cens_ij_env) - 1), N = n_samp
          )
        }
        y_pred_VT[mask_ij_env][mask_inner_ij_env] <- rowMeans(samp_ij_env_VT)[mask_inner_ij_env]
        sd_pred_VT[mask_ij_env][mask_inner_ij_env] <-
          apply(samp_ij_env_VT, 1, sd)[mask_inner_ij_env] / sqrt(n_samp)
        y_samp_VT[mask_ij_env, ][mask_inner_ij_env, ] <- samp_ij_env_VT[mask_inner_ij_env, ]
        end_time <- Sys.time()
        time_VT <- time_VT + difftime(end_time, bgn_time, units = "secs")[[1]]
      }

      if (run_TN) {
        bgn_time <- Sys.time()
        set.seed(123)
        cat("TN sampling at i =", i, "j =", j, "\n")
        covmat_ij_env <- getFromNamespace(cov_name, "GpGp")(cov_parms,
          locs_ij_env)
        if (sum(mask_cens_ij_env) == sum(mask_ij_env)) {
          cond_mean_ij_env_cens <- rep(0, sum(mask_ij_env))
          cond_covmat_ij_env_cens <- covmat_ij_env
        } else {
          cond_mean_ij_env_cens <- as.vector(covmat_ij_env[mask_cens_ij_env, !mask_cens_ij_env] %*%
            solve(
              covmat_ij_env[!mask_cens_ij_env, !mask_cens_ij_env],
              y_obs_ij_env[!mask_cens_ij_env]
            ))
          cond_covmat_ij_env_cens <- covmat_ij_env[mask_cens_ij_env, mask_cens_ij_env] -
            covmat_ij_env[mask_cens_ij_env, !mask_cens_ij_env] %*%
            solve(covmat_ij_env[!mask_cens_ij_env, !mask_cens_ij_env]) %*%
            covmat_ij_env[!mask_cens_ij_env, mask_cens_ij_env]
        }
        cond_covmat_ij_env_cens[lower.tri(cond_covmat_ij_env_cens)] <-
          t(cond_covmat_ij_env_cens)[lower.tri(cond_covmat_ij_env_cens)]
        samp_ij_env_TN <- t(TruncatedNormal::rtmvnorm(
          n_samp, cond_mean_ij_env_cens,
          cond_covmat_ij_env_cens, rep(-Inf, sum(mask_cens_ij_env)),
          levl_cens_ij_env[mask_cens_ij_env]
        ))
        locs_cens_ij_env <- locs_ij_env[mask_cens_ij_env, , drop = FALSE]
        mask_cens_inner_ij_env <- locs_cens_ij_env[, 1] >= (i - 1) * 0.2 &
          locs_cens_ij_env[, 1] <= i * 0.2 &
          locs_cens_ij_env[, 2] >= (j - 1) * 0.2 & locs_cens_ij_env[, 2] <= j * 0.2
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
      file = paste0("results/sim_data_VT_scene", scene_ID, ".RData")
    )
    cat(
      "RMSE for VT is ",
      sqrt(mean((y[mask_cens] - y_pred_VT[mask_cens])^2)), "\n"
    )
    cat(
      "NLL for VT is ",
      -mean(dnorm(y[mask_cens],
        mean = y_pred_VT[mask_cens],
        sd = sd_pred_VT[mask_cens]
      )), "\n"
    )
    cat(
      "CRPS for VT is ",
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
      file = paste0("results/sim_data_TN_scene", scene_ID, ".RData")
    )
    cat(
      "RMSE for TN is ",
      sqrt(mean((y[mask_cens] - y_pred_TN[mask_cens])^2)), "\n"
    )
    cat(
      "NLL for TN is ",
      -mean(dnorm(y[mask_cens],
        mean = y_pred_TN[mask_cens],
        sd = sd_pred_TN[mask_cens]
      )), "\n"
    )
    cat(
      "CRPS for TN is ",
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
    scene_ID, "_m", m, "_order", reorder, ".RData"
  ))
  load(paste0("results/sim_data_VT_scene", scene_ID, ".RData"))
  load(paste0("results/sim_data_TN_scene", scene_ID, ".RData"))
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
