library(GpGp)
library(VeccTMVN)

rm(list = ls())
scene_ID <- 1
m <- 30
n_samp <- 50
subset_size <- 2500
n_burn <- 20000
n_iter_MC <- 25000
thin <- 5

args <- commandArgs(trailingOnly = TRUE)
k <- if (length(args)) as.integer(args[1]) else 1

if (!requireNamespace("CensSpBayes", quietly = TRUE)) {
  stop("This script requires the CensSpBayes package.")
}

# Common simulated data ---------------------------------------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs_SNN <- y
y_obs_SNN[mask_cens] <- NA
y_obs_CSB <- y
y_obs_CSB[mask_cens] <- cens_ub[mask_cens]
source("../utils/score_output.R")

# CSB ---------------------------------------------------------------------
set.seed(123)
bgn_time <- Sys.time()
inla_mats <- CensSpBayes::create_inla_mats(
  S = locs, S.pred = locs[mask_cens, , drop = FALSE],
  offset = c(0.01, 0.2), cutoff = 0.05, max.edge = c(0.01, 0.1)
)
ret_CSB <- CensSpBayes::CensSpBayes(
  Y = y_obs_CSB, S = locs, X = matrix(1, n, 1),
  cutoff.Y = cens_ub,
  S.pred = locs[mask_cens, , drop = FALSE],
  X.pred = matrix(1, sum(mask_cens), 1),
  inla.mats = inla_mats, rho.init = 0.1, rho.upper = 5,
  iters = n_iter_MC, burn = n_burn, thin = thin, ret_samp = TRUE
)
comp_time_CSB <- difftime(Sys.time(), bgn_time, units = "secs")[[1]]
y_cens_samp_CSB <- tail(ret_CSB$Y.pred.samp, sum(mask_cens))
score_output(
  y_cens_samp_CSB,
  y[mask_cens],
  comp_time_CSB,
  scene_ID = scene_ID,
  method = "CSB",
  parms = "unknown"
)

rm(inla_mats, ret_CSB, y_cens_samp_CSB, covmat)
gc()

# SNN with fitted covariance and maximin ordering -------------------------
order_est <- head(GpGp::order_maxmin(locs), min(subset_size, n))
nll <- function(log_parms) {
  parms <- exp(log_parms)
  set.seed(123)
  -VeccTMVN::loglk_censor_MVN(
    locs[order_est, , drop = FALSE], which(mask_cens[order_est]),
    y_obs_SNN[order_est], cens_ub[order_est],
    covName = "matern15_isotropic", covParms = parms, m = m
  )
}

bgn_time <- Sys.time()
cov_parms_est <- exp(optim(log(c(0.5, 0.01, 0.01)), nll,
  control = list(maxit = 200)
)$par)

ord <- GpGp::order_maxmin(locs)
K_order <- GpGp::matern15_isotropic(cov_parms_est, locs[ord, , drop = FALSE])
y_samp_order <- lapply(seq_len(n_samp), function(seed) {
  nntmvn::rptmvn(
    y_obs_SNN[ord], cens_lb[ord], cens_ub[ord], mask_cens[ord],
    m = m, covmat = K_order, locs = locs[ord, , drop = FALSE],
    ordering = 0, seed = seed
  )
})
comp_time_SNN <- difftime(Sys.time(), bgn_time, units = "secs")[[1]]
reverse_order <- order(ord)
y_samp_SNN <- matrix(unlist(y_samp_order), n, n_samp)[reverse_order, , drop = FALSE]
rm(K_order, y_samp_order)
gc()

score_output(
  y_samp_SNN[mask_cens],
  y[mask_cens],
  comp_time_SNN,
  scene_ID = scene_ID,
  m = m,
  method = "SNN_order_maximin",
  parms = "unknown"
)
