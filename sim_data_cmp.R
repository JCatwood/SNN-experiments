library(mvtnorm)
library(RANN)
library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)

# input args ---------------
rm(list = ls())
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  scene_ID <- as.integer(args[1])
  m <- as.integer(args[2])
  reorder <- as.integer((args[3]))
  run_seq_Vecc <- args[4] == "T" || args[4] == "True" || args[4] == "true" ||
    args[4] == "TRUE"
  run_CB <- args[5] == "T" || args[5] == "True" || args[5] == "true" ||
    args[5] == "TRUE"
  run_est_seq_Vecc <- args[6] == "T" || args[6] == "True" || args[6] == "true" ||
    args[6] == "TRUE"
  subset_size <- as.integer((args[7]))
  plot_heatmap <- args[8] == "T" || args[8] == "True" || args[8] == "true" ||
    args[8] == "TRUE"
  run_VT <- args[9] == "T" || args[9] == "True" || args[9] == "true" ||
    args[9] == "TRUE"
  run_TN <- args[10] == "T" || args[10] == "True" || args[10] == "true" ||
    args[10] == "TRUE"
} else {
  scene_ID <- 1
  m <- 30
  reorder <- 0 # 0 no reorder, 1 maximin
  run_seq_Vecc <- FALSE
  run_CB <- FALSE
  run_est_seq_Vecc <- FALSE
  subset_size <- 2500
  plot_heatmap <- FALSE
  run_VT <- TRUE
  run_TN <- TRUE
}

# sim settings ------------------------
n_samp <- 50
n_burn <- 20000
n_iter_MC <- 25000
thin <- 5
if (scene_ID == 1) {
  tmp_vec <- seq(from = 0, to = 1, length.out = 100)
  locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
  cov_func <- GpGp::matern15_isotropic
  cov_parms <- c(1.0, 0.03, 0.0001)
  cov_name <- "matern15_isotropic"
  covmat <- cov_func(cov_parms, locs)
  cat("Generating GP ...", "\n")
  y <- as.vector(mvtnorm::rmvnorm(1, sigma = covmat))
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
  library(nntmvn)
  ncores <- 4
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  bgn_time <- Sys.time()
  if (reorder == 0) {
    order <- 1:n
    y_obs_order <- y_obs
    levl_cens_order <- levl_cens
    mask_cens_order <- mask_cens
    NN_order <- RANN::nn2(locs, k = m + 1)[[1]]
    covmat_order <- covmat
  } else if (reorder == 1) {
    order <- GpGp::order_maxmin(locs)
    y_obs_order <- y_obs[order]
    levl_cens_order <- levl_cens[order]
    mask_cens_order <- mask_cens[order]
    NN_order <- RANN::nn2(locs[order, , drop = FALSE], k = m + 1)[[1]]
    covmat_order <- covmat[order, order]
  }
  y_samp_seq_Vecc_order <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
    nntmvn::seq_Vecc_samp_func(y_obs_order, rep(-Inf, n), levl_cens_order,
      mask_cens_order, NN_order,
      covmat = covmat_order,
      seed = i
    )
  }
  end_time <- Sys.time()
  stopCluster(cl)
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
    "NLL for est seq Vecc is ",
    -mean(dnorm(y[mask_cens],
      mean = y_pred_cens_seq_Vecc,
      sd = sd_cens_seq_Vecc
    )), "\n"
  )
}

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
  y_samp_CB <- CensSpBayes::CensSpBayes(
    Y = y_obs, S = locs, X = X.obs,
    cutoff.Y = levl_cens,
    S.pred = locs[mask_cens, ], X.pred = X.pred,
    inla.mats = inla.mats,
    rho.init = 0.1, rho.upper = 5,
    iters = n_iter_MC, burn = n_burn, thin = thin, ret_samp = TRUE
  )
  end_time <- Sys.time()
  time_CB <- difftime(end_time, bgn_time, units = "secs")[[1]]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(time_CB, y_samp_CB, file = paste0(
    "results/sim_data_CB_scene",
    scene_ID, ".RData"
  ))
  cat(
    "RMSE for CB is ",
    sqrt(mean((y[mask_cens] - y_samp_CB$Y.pred.posmean)^2)),
    "\n"
  )
  cat(
    "NLL for CB is ",
    -mean(dnorm(y[mask_cens],
      mean = y_samp_CB$Y.pred.posmean,
      sd = sqrt(y_samp_CB$Y.pred.posvar /
        floor((n_iter_MC - n_burn) / thin))
    )),
    "\n"
  )
}

# nntmvn with est. parms -------------------------------
if (run_est_seq_Vecc) {
  library(VeccTMVN)
  library(nntmvn)
  if (reorder == 0) {
    if (subset_size < n) {
      order <- floor(seq(from = 1, to = n, length.out = subset_size))
    } else {
      order <- 1:n
    }
  } else if (reorder == 1) {
    order <- GpGp::order_maxmin(locs)
    if (subset_size < n) {
      order <- order[1:subset_size]
    }
  }
  y_obs_order <- y_obs[order]
  levl_cens_order <- levl_cens[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  NN_order <- RANN::nn2(locs_order, k = m + 1)[[1]]
  if (sum(mask_cens_order) == min(subset_size, n)) {
    stop("No observed responses in this subset\n")
  }

  nll_func <- function(log_parms) {
    parms <- exp(log_parms)
    set.seed(123)
    ll <- VeccTMVN::loglk_censor_MVN(locs_order, which(mask_cens_order),
      y_obs_order, levl_cens_order,
      covName = cov_name, covParms = parms, m = m
    )
    cat("Parms = ", parms, ", nll = ", -ll, "\n")
    -ll
  }
  logparms_init <- log(c(0.5, 0.01, 0.01))
  bgn_time <- Sys.time()
  cov_parms_est <- exp(optim(logparms_init, nll_func,
    control = list(maxit = 200)
  )$par)
  end_time <- Sys.time()
  time_parm_est_seq_Vecc <- difftime(end_time, bgn_time, units = "secs")[[1]]
  cat("Parameter estimation used ", time_parm_est_seq_Vecc, " seconds\n")
  cat("Estimated parameters are ", cov_parms_est, "\n")

  ncores <- 4
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  bgn_time <- Sys.time()
  if (reorder == 0) {
    order <- 1:n
  } else if (reorder == 1) {
    order <- GpGp::order_maxmin(locs)
  }
  y_obs_order <- y_obs[order]
  levl_cens_order <- levl_cens[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  NN_order <- RANN::nn2(locs_order, k = m + 1)[[1]]
  covmat_order <- getFromNamespace(cov_name, "GpGp")(cov_parms_est, locs_order)
  y_samp_est_seq_Vecc_order <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
    nntmvn::seq_Vecc_samp_func(y_obs_order, rep(-Inf, n), levl_cens_order,
      mask_cens_order, NN_order,
      covmat = covmat_order,
      seed = i
    )
  }
  end_time <- Sys.time()
  stopCluster(cl)
  time_est_seq_Vecc <- difftime(end_time, bgn_time, units = "secs")[[1]]
  rev_order <- 1:n
  rev_order[order] <- 1:n
  y_samp_est_seq_Vecc <- matrix(unlist(y_samp_est_seq_Vecc_order), n,
    n_samp,
    byrow = FALSE
  )[rev_order, , drop = FALSE]
  y_pred_est_seq_Vecc <- rowMeans(y_samp_est_seq_Vecc)
  y_pred_cens_est_seq_Vecc <- y_pred_est_seq_Vecc[mask_cens]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(time_parm_est_seq_Vecc, time_est_seq_Vecc, y_samp_est_seq_Vecc,
    file = paste0(
      "results/sim_data_est_seq_Vecc_scene",
      scene_ID, "_m", m, "_order", reorder, "_subset", subset_size, ".RData"
    )
  )
  cat(
    "RMSE for est seq Vecc is ",
    sqrt(mean((y[mask_cens] - y_pred_cens_est_seq_Vecc)^2)), "\n"
  )
  sd_cens_est_seq_Vecc <- apply(y_samp_est_seq_Vecc, 1, sd)[mask_cens] /
    sqrt(n_samp)
  cat(
    "NLL for est seq Vecc is ",
    -mean(dnorm(y[mask_cens],
      mean = y_pred_cens_est_seq_Vecc,
      sd = sd_cens_est_seq_Vecc
    )), "\n"
  )
}

# VeccTMVN and TruncatedNormal -----------------------
if (run_VT || run_TN) {
  if (run_VT) {
    y_pred_VT <- y_obs
    sd_pred_VT <- rep(0, length(y_obs))
    y_samp_VT <- y_obs
    time_VT <- 0
  }
  if (run_TN) {
    y_pred_TN <- y_obs
    sd_pred_TN <- rep(0, length(y_obs))
    y_samp_TN <- y_obs
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
        y_samp_VT[mask_ij_env][mask_inner_ij_env] <- samp_ij_env_VT[mask_inner_ij_env, 1]
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
        y_samp_TN[mask_ij_env][mask_inner_ij_env & mask_cens_ij_env] <-
          samp_ij_env_TN[mask_cens_inner_ij_env, 1]
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
  }
}

# heatmap -------------------------
if (plot_heatmap) {
  library(fields)
  load(paste0(
    "results/sim_data_CB_scene",
    scene_ID, ".RData"
  ))
  load(paste0(
    "results/sim_data_seq_Vecc_scene",
    scene_ID, "_m", m, "_order", reorder, ".RData"
  ))
  load(paste0("results/sim_data_VT_scene", scene_ID, ".RData"))
  load(paste0("results/sim_data_TN_scene", scene_ID, ".RData"))
  y_samp_CB_tmp <- y_obs
  y_samp_CB_tmp[mask_cens] <- y_samp_CB$Y.pred.samp[, ncol(y_samp_CB$Y.pred.samp)]
  y_samp_CB <- y_samp_CB_tmp
  y_samp_seq_Vecc <- y_samp_seq_Vecc[, 1]
  rm(y_samp_CB_tmp)
  zlim <- range(y_samp_CB, y_samp_seq_Vecc, y_samp_TN, y_samp_VT)
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
    file = paste0("plots/sim_data_CB_scene", scene_ID, ".pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_CB, sqrt(n), sqrt(n)), zlim = zlim)
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
