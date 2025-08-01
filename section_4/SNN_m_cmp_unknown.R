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
reorder <- 1 # 0 no reorder, 1 maximin
n_samp <- 50 # samples generated for posterior inference
subset_size <- 2500
run_est_SNN <- TRUE
args <- commandArgs(trailingOnly = TRUE)
# use command line args when running in batch on clusters
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2]) # simulation scenario ID
  reorder <- as.integer(args[3]) # 0 no reorder, 1 maximin
}

# data simulation ----------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- NA
if (!exists("cov_name")) {
  cov_name <- "matern15_isotropic"
}

source("../utils/score_output.R")

# nntmvn ---------------------------------------
for (m in m_seq) {
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
  cens_ub_order <- cens_ub[order]
  cens_lb_order <- cens_lb[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  if (sum(mask_cens_order) == min(subset_size, n)) {
    stop("No observed responses in this subset\n")
  }
  nll_func <- function(log_parms) {
    parms <- exp(log_parms)
    set.seed(123)
    ll <- VeccTMVN::loglk_censor_MVN(locs_order, which(mask_cens_order),
      y_obs_order, cens_ub_order,
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
  time_parm_est_SNN <- difftime(end_time, bgn_time, units = "secs")[[1]]
  cat("Parameter estimation used ", time_parm_est_SNN, " seconds\n")
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
  cens_ub_order <- cens_ub[order]
  cens_lb_order <- cens_lb[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  covmat_order <- getFromNamespace(cov_name, "GpGp")(cov_parms_est, locs_order)
  bgn_time <- Sys.time()
  y_samp_est_SNN_order <- lapply(1:n_samp, function(seed_id) {
    nntmvn::rptmvn(y_obs_order, cens_lb_order, cens_ub_order,
      mask_cens_order,
      m = m,
      covmat = covmat_order, locs = locs_order,
      seed = seed_id
    )
  })
  end_time <- Sys.time()
  cat("SNN sampling is done\n")
  time_est_SNN <- difftime(end_time, bgn_time, units = "secs")[[1]]
  rev_order <- 1:n
  rev_order[order] <- 1:n
  y_samp_est_SNN <- matrix(unlist(y_samp_est_SNN_order), n,
    n_samp,
    byrow = FALSE
  )[rev_order, , drop = FALSE]

  if (reorder == 0) {
    method <- "SNN"
  } else if (reorder == 1) {
    method <- "SNN_order_maximin"
  } else {
    stop("Unknown reorder\n")
  }

  score_output(y_samp_est_SNN[mask_cens, ], y[mask_cens],
    time_est_SNN + time_parm_est_SNN,
    scene_ID = scene_ID, method = method, parms = "unknown"
  )
}
