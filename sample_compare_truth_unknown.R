library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(scoringRules)

# simulation settings ---------------
rm(list = ls())
set.seed(123)
# general
scene_ID <- 1 # Matern 1.5 kernel
m <- 30 # number of nearest neighbors
reorder <- 0 # 0 no reorder, 1 maximin
n_samp <- 50 # samples generated for posterior inference
subset_size <- 2500
run_est_seq_Vecc <- TRUE
run_CB <- FALSE
use_parallel <- FALSE
plot_heatmap <- FALSE
# CensSpBayes
n_burn <- 20000
n_iter_MC <- 25000
thin <- 5

# data simulation ----------------------
if (scene_ID == 1) {
  tmp_vec <- seq(from = 0, to = 1, length.out = 50)
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
if (scene_ID == 2) {
  tmp_vec <- seq(from = 0, to = 1, length.out = 100)
  locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
}

mask_cens <- y < levl_cens
y_obs <- y
y_obs[mask_cens] <- levl_cens[mask_cens]

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
  cat(
    "CRPS for CB is ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens],
      dat = y_samp_CB$Y.pred.samp
    )), "\n"
  )
}

# nntmvn with est. parms -------------------------------
if (run_est_seq_Vecc) {
  ## covariance parameter estimation --------------------------------
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
  ## sampling ----------------------
  if (reorder == 0) {
    order <- 1:n
  } else if (reorder == 1) {
    order <- GpGp::order_maxmin(locs)
  }
  cat("SNN sampling...\n")
  bgn_time <- Sys.time()
  y_obs_order <- y_obs[order]
  levl_cens_order <- levl_cens[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  covmat_order <- getFromNamespace(cov_name, "GpGp")(cov_parms_est, locs_order)
  if (use_parallel) {
    ncores <- 4
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    if (reorder == 0) {
      order <- 1:n
    } else if (reorder == 1) {
      order <- GpGp::order_maxmin(locs)
    }
    y_samp_est_seq_Vecc_order <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
      nntmvn::rtmvn_snn(y_obs_order, rep(-Inf, n), levl_cens_order,
        mask_cens_order,
        m = m,
        covmat = covmat_order, locs = locs_order,
        seed = i
      )
    }
    stopCluster(cl)
  } else {
    y_samp_est_seq_Vecc_order <- lapply(1:n_samp, function(seed_id) {
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
  cat(
    "CRPS for est seq Vecc is ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens],
      dat = y_samp_est_seq_Vecc[mask_cens, , drop = FALSE]
    )), "\n"
  )
}
# heatmap -------------------------
if (plot_heatmap) {
  library(fields)
  load(paste0(
    "results/sim_data_CB_scene",
    scene_ID, ".RData"
  ))
  load(paste0(
    "results/sim_data_est_seq_Vecc_scene",
    scene_ID, "_m", m, "_order", reorder, ".RData"
  ))
  y_samp_CB_tmp <- y_obs
  y_samp_CB_tmp[mask_cens] <- y_samp_CB$Y.pred.samp[, ncol(y_samp_CB$Y.pred.samp)]
  y_samp_CB <- y_samp_CB_tmp
  y_samp_est_seq_Vecc <- y_samp_est_seq_Vecc[, 1]
  rm(y_samp_CB_tmp)
  zlim <- range(y_samp_CB, y_samp_est_seq_Vecc)
  if (!file.exists("plots")) {
    dir.create("plots")
  }
  pdf(
    file = paste0("plots/sim_heatmap_truth_unknown_", scene_ID, "_true.pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y, sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0(
      "plots/sim_data_est_seq_Vecc_scene", scene_ID, "_m", m,
      "_order", reorder, ".pdf"
    ),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_est_seq_Vecc, sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0("plots/sim_data_CB_scene", scene_ID, ".pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_CB, sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
}
