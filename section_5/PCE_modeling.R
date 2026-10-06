library(VeccTMVN)
library(sf)
library(spData)
library(GpGp)

rm(list = ls())
set.seed(123)

run_parm_est <- FALSE
run_VeccTMVN <- TRUE
run_TN <- TRUE
run_CB <- FALSE
run_seq_Vecc <- TRUE
run_seq_Vecc_all <- TRUE
run_CB_all <- FALSE
n_national <- 3 # independent SNN draws; retain the same number for CSB

# CensSpBayes is optional and is checked only when a CSB experiment is enabled.
if ((run_CB || run_CB_all) &&
    !requireNamespace("CensSpBayes", quietly = TRUE)) {
  stop(
    "run_CB or run_CB_all requires the optional CensSpBayes package.",
    call. = FALSE
  )
}

load("PCE.RData")
# extract raw data -------------------------------
y <- log(data.PCE.censored$result_va + 1e-8)
b_censor <- log(data.PCE.censored$detection_level + 1e-8)
mask_censor <- data.PCE.censored$left_censored
locs <- cbind(
  data.PCE.censored$lon, data.PCE.censored$lat,
  data.PCE.censored$startDate
)
n <- nrow(locs)
d <- ncol(locs)
locs_names <- c("lon", "lat", "date")
colnames(locs) <- locs_names
# standardize -------------------------------
b_scaled <- (b_censor - mean(y, na.rm = T)) / sd(y, na.rm = T)
y_scaled <- (y - mean(y, na.rm = T)) / sd(y, na.rm = T)
locs_scaled <- locs
for (j in 1:d) {
  locs_scaled[, j] <- (locs_scaled[, j] - min(locs[, j])) /
    (max(locs[, j]) - min(locs[, j]))
}
# ordering and subsetting --------------------
subset_size <- round(n * 0.2)
order <- GpGp::order_maxmin(locs[, 1:2])[1:subset_size]
# model fitting -------------------------------
cov_name <- "matern15_scaledim"
if (run_parm_est) {
  logcovparms_init <- log(c(1, 0.01, 0.01, 0.01, 0.005)) # var, range1, range2, range3, nugget
  neglk_func <- function(logcovparms, ...) {
    covparms <- exp(logcovparms)
    set.seed(123)
    negloglk <- -loglk_censor_MVN(
      locs_scaled[order, , drop = FALSE], which(mask_censor[order]),
      y_scaled[order], b_scaled[order], cov_name,
      covparms, ...
    )
    cat("covparms is", covparms, "\n")
    cat("Neg loglk is", negloglk, "\n")
    return(negloglk)
  }
  opt_obj <- optim(
    par = logcovparms_init, fn = neglk_func,
    control = list(trace = 1, maxit = 200), m = 50, NLevel2 = 1e3
  )
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(opt_obj, file = paste0(
    "results/PCE_modeling.RData"
  ))
}
# Find State information for locs -----------------------------------
lonlat_to_state <- function(locs) {
  ## State DF
  states <- spData::us_states
  ## Convert points data.frame to an sf POINTS object
  pts <- st_as_sf(locs, coords = 1:2, crs = 4326)
  ## Transform spatial data to some planar coordinate system
  ## (e.g. Web Mercator) as required for geometric operations
  states <- st_transform(states, crs = 3857)
  pts <- st_transform(pts, crs = 3857)
  ## Find names of state (if any) intersected by each point
  state_names <- states[["NAME"]]
  ii <- vapply(st_intersects(pts, states), function(i) if (length(i)) i[1] else NA_integer_, integer(1))
  state_names[ii]
}
state_names <- lonlat_to_state(data.frame(locs))
ind_Texas <- which(state_names == "Texas")
# Define a enveloping box for Texas -------------------------------------------
## TX boundary lat: 25.83333 to 36.5 lon: -93.51667 to -106.6333
ind_Texas_box <- which((locs[, 1] < -93.02) & (locs[, 1] > -107.13) &
  (locs[, 2] > 25.33) & (locs[, 2] < 37))
# Sample with VeccTMVN over TX ------------------
if (run_VeccTMVN) {
  load("results/PCE_modeling.RData")
  ind_Texas_big <- ind_Texas
  ind_obs <- which(!is.na(y))
  ind_Texas_big <- union(ind_Texas_big, ind_obs)
  ind_censor_Texas_big <- which(is.na(y[ind_Texas_big]))
  locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
  y_scaled_Texas_big <- y_scaled[ind_Texas_big]
  b_scaled_Texas_big <- b_scaled[ind_Texas_big]
  ind_obs_tmp <- which(!is.na(y_scaled_Texas_big))
  n_samp <- 3
  n_cens <- sum(is.na(y_scaled_Texas_big))
  covparms <- exp(opt_obj$par)
  m <- 50
  time_bgn <- Sys.time()
  covmat_Texas_big <- getFromNamespace(cov_name, "GpGp")(covparms,
    locs_scaled_Texas_big)
  cond_mean_TX_cens <- as.vector(covmat_Texas_big[-ind_obs_tmp, ind_obs_tmp] %*%
    solve(
      covmat_Texas_big[ind_obs_tmp, ind_obs_tmp],
      y_scaled_Texas_big[ind_obs_tmp]
    ))
  cond_covmat_TX_cens <- covmat_Texas_big[-ind_obs_tmp, -ind_obs_tmp] -
    covmat_Texas_big[-ind_obs_tmp, ind_obs_tmp] %*%
    solve(covmat_Texas_big[ind_obs_tmp, ind_obs_tmp]) %*%
    covmat_Texas_big[ind_obs_tmp, -ind_obs_tmp]
  cond_covmat_TX_cens[lower.tri(cond_covmat_TX_cens)] <-
    t(cond_covmat_TX_cens)[lower.tri(cond_covmat_TX_cens)]
  samp_TX_VT <- outer(y_scaled_Texas_big, rep(1, n_samp))
  time_TX_VT <- system.time(
    samp_TX_VT[-ind_obs_tmp, ] <- VeccTMVN::mvrandn(
      lower = rep(-Inf, n_cens), upper = b_scaled_Texas_big[-ind_obs_tmp],
      mean = cond_mean_TX_cens, sigma = cond_covmat_TX_cens, N = n_samp, m = m
    )
  )[[3]]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(samp_TX_VT, time_TX_VT, covparms, cov_name, ind_Texas_big,
    file = "results/PCE_samp_VT.RData"
  )
}
# SNN method for TX --------------------
if (run_seq_Vecc) {
  library(nntmvn)
  library(doParallel)
  load("results/PCE_modeling.RData")
  ind_Texas_big <- ind_Texas
  ind_obs <- which(!is.na(y))
  ind_Texas_big <- union(ind_Texas_big, ind_obs)
  ind_censor_Texas_big <- which(is.na(y[ind_Texas_big]))
  locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
  y_scaled_Texas_big <- y_scaled[ind_Texas_big]
  b_scaled_Texas_big <- b_scaled[ind_Texas_big]
  ind_obs_tmp <- which(!is.na(y_scaled_Texas_big))
  n_samp <- 3
  n_cens <- sum(is.na(y_scaled_Texas_big))
  covparms <- exp(opt_obj$par)
  m <- 50
  time_bgn <- Sys.time()
  covmat_Texas_big <- getFromNamespace(cov_name, "GpGp")(covparms,
    locs_scaled_Texas_big)
  cond_mean_TX_cens <- as.vector(covmat_Texas_big[-ind_obs_tmp, ind_obs_tmp] %*%
    solve(
      covmat_Texas_big[ind_obs_tmp, ind_obs_tmp],
      y_scaled_Texas_big[ind_obs_tmp]
    ))
  cond_covmat_TX_cens <- covmat_Texas_big[-ind_obs_tmp, -ind_obs_tmp] -
    covmat_Texas_big[-ind_obs_tmp, ind_obs_tmp] %*%
    solve(covmat_Texas_big[ind_obs_tmp, ind_obs_tmp]) %*%
    covmat_Texas_big[ind_obs_tmp, -ind_obs_tmp]
  cond_covmat_TX_cens[lower.tri(cond_covmat_TX_cens)] <-
    t(cond_covmat_TX_cens)[lower.tri(cond_covmat_TX_cens)]
  samp_seq_Vecc <- outer(y_scaled_Texas_big, rep(1, n_samp))
  ncores <- 4
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  samp_seq_Vecc_sub <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
    nntmvn::rtmvn(
      cens_lb = rep(-Inf, n_cens) - cond_mean_TX_cens,
      cens_ub = b_scaled_Texas_big[-ind_obs_tmp] - cond_mean_TX_cens,
      m = m, covmat = cond_covmat_TX_cens, ordering = 2,
      locs = locs_scaled_Texas_big[-ind_obs_tmp, ], seed = i
    ) + cond_mean_TX_cens
  }
  stopCluster(cl)
  time_end <- Sys.time()
  time_seq_Vecc <- difftime(time_end, time_bgn, units = "secs")[[1]]
  samp_seq_Vecc_sub <- matrix(unlist(samp_seq_Vecc_sub),
    length(samp_seq_Vecc_sub[[1]]),
    n_samp,
    byrow = FALSE
  )
  samp_seq_Vecc[-ind_obs_tmp, ] <- samp_seq_Vecc_sub
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(samp_seq_Vecc, time_seq_Vecc, covparms, cov_name, ind_Texas_big,
    file = "results/PCE_samp_seq_Vecc.RData"
  )
}
# SNN method for all locations --------------------
if (run_seq_Vecc_all) {
  library(nntmvn)
  load("results/PCE_modeling.RData")
  covparms <- exp(opt_obj$par)
  m <- 50
  locs_scaled_twice <- locs_scaled %*% diag(1 / covparms[2:4])
  covparms_tmp <- covparms
  covparms_tmp[2:4] <- 1
  time_bgn <- Sys.time()
  samp_seq_Vecc_all <- vapply(seq_len(n_national), function(s) nntmvn::rptmvn(
    y_scaled, rep(-Inf, n), b_scaled, is.na(y_scaled), m,
    locs = locs_scaled_twice, cov_name = cov_name,
    cov_parm = covparms_tmp, seed = s, ordering = 2
  ), numeric(n))
  time_end <- Sys.time()
  time_seq_Vecc_all <- difftime(time_end, time_bgn, units = "secs")[[1]]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(samp_seq_Vecc_all, time_seq_Vecc_all, covparms, cov_name,
    file = "results/PCE_samp_seq_Vecc_all.RData"
  )
}
# CensSpBayes method --------------------
if (run_CB) {
  ind_Texas_big <- ind_Texas
  ind_obs <- which(!is.na(y))
  ind_Texas_big <- union(ind_Texas_big, ind_obs)
  locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
  y_scaled_Texas_big <- y_scaled[ind_Texas_big]
  b_scaled_Texas_big <- b_scaled[ind_Texas_big]
  mask_cens_Texas_big <- is.na(y_scaled_Texas_big)
  y_obs_scaled_Texas_big <- y_scaled_Texas_big
  y_obs_scaled_Texas_big[mask_cens_Texas_big] <-
    b_scaled_Texas_big[mask_cens_Texas_big]
  n_burn <- 20000
  n_iter_MC <- 20015
  thin <- 5
  ## Sample at locations given by `ind_Texas_big` using CB -----------------
  time_bgn <- Sys.time()
  inla.mats <- CensSpBayes::create_inla_mats(
    S = locs_scaled_Texas_big[, 1:2], # 3D mesh produced error
    S.pred = locs_scaled_Texas_big[mask_cens_Texas_big, 1:2],
    offset = c(0.01, 0.2),
    cutoff = 0.05,
    max.edge = c(0.01, 0.1)
  )
  X.obs <- matrix(1, nrow(locs_scaled_Texas_big), 1)
  X.pred <- matrix(1, sum(mask_cens_Texas_big), 1)
  y_samp_CB <- CensSpBayes::CensSpBayes(
    Y = y_obs_scaled_Texas_big, S = locs_scaled_Texas_big[, 1:2], X = X.obs,
    cutoff.Y = b_scaled_Texas_big,
    S.pred = locs_scaled_Texas_big[mask_cens_Texas_big, 1:2], X.pred = X.pred,
    inla.mats = inla.mats,
    rho.init = 0.1, rho.upper = 5,
    iters = n_iter_MC, burn = n_burn, thin = thin, ret_samp = TRUE
  )
  time_end <- Sys.time()
  time_TX_CB <- difftime(time_end, time_bgn, units = "secs")[[1]]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(y_samp_CB, time_TX_CB, ind_Texas_big,
    file = "results/PCE_samp_CB.RData"
  )
}
if (run_CB_all) {
  y_obs_scaled <- y_scaled
  y_obs_scaled[mask_censor] <- b_scaled[mask_censor]
  n_burn <- 20000
  n_iter_MC <- n_burn + 5 * n_national
  thin <- 5
  ## Sample at all locations using CB -----------------
  time_bgn <- Sys.time()
  inla.mats <- CensSpBayes::create_inla_mats(
    S = locs_scaled[, 1:2],
    S.pred = locs_scaled[mask_censor, 1:2],
    offset = c(0.01, 0.2),
    cutoff = 0.05,
    max.edge = c(0.01, 0.1)
  )
  X.obs <- matrix(1, nrow(locs_scaled), 1)
  X.pred <- matrix(1, sum(mask_censor), 1)
  y_samp_CB_all <- CensSpBayes::CensSpBayes(
    Y = y_obs_scaled, S = locs_scaled[, 1:2], X = X.obs,
    cutoff.Y = b_scaled,
    S.pred = locs_scaled[mask_censor, 1:2], X.pred = X.pred,
    inla.mats = inla.mats,
    rho.init = 0.1, rho.upper = 5,
    iters = n_iter_MC, burn = n_burn, thin = thin, ret_samp = TRUE
  )
  time_end <- Sys.time()
  time_TX_CB_all <- difftime(time_end, time_bgn, units = "secs")[[1]]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(y_samp_CB_all, time_TX_CB_all, file = "results/PCE_samp_CB_all.RData")
}
# TruncatedNormal method --------------------
if (run_TN) {
  library(TruncatedNormal)
  load("results/PCE_modeling.RData")
  ind_Texas_big <- ind_Texas
  ind_obs <- which(!is.na(y))
  ind_Texas_big <- union(ind_Texas_big, ind_obs)
  locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
  y_scaled_Texas_big <- y_scaled[ind_Texas_big]
  b_scaled_Texas_big <- b_scaled[ind_Texas_big]
  ind_obs_tmp <- which(!is.na(y_scaled_Texas_big))
  n_samp <- 3
  n_cens <- sum(is.na(y_scaled_Texas_big))
  covparms <- exp(opt_obj$par)
  ## Sample at locations given by `ind_Texas_big` using TN -----------------
  time_bgn <- Sys.time()
  covmat_Texas_big <- getFromNamespace(cov_name, "GpGp")(covparms,
    locs_scaled_Texas_big)
  cond_mean_TX_cens <- as.vector(covmat_Texas_big[-ind_obs_tmp, ind_obs_tmp] %*%
    solve(
      covmat_Texas_big[ind_obs_tmp, ind_obs_tmp],
      y_scaled_Texas_big[ind_obs_tmp]
    ))
  cond_covmat_TX_cens <- covmat_Texas_big[-ind_obs_tmp, -ind_obs_tmp] -
    covmat_Texas_big[-ind_obs_tmp, ind_obs_tmp] %*%
    solve(covmat_Texas_big[ind_obs_tmp, ind_obs_tmp]) %*%
    covmat_Texas_big[ind_obs_tmp, -ind_obs_tmp]
  cond_covmat_TX_cens[lower.tri(cond_covmat_TX_cens)] <-
    t(cond_covmat_TX_cens)[lower.tri(cond_covmat_TX_cens)]
  samp_TX_TN <- TruncatedNormal::rtmvnorm(
    n_samp, cond_mean_TX_cens,
    cond_covmat_TX_cens, rep(-Inf, n_cens), b_scaled_Texas_big[-ind_obs_tmp]
  )
  time_end <- Sys.time()
  time_TX_TN <- difftime(time_end, time_bgn, units = "secs")[[1]]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(samp_TX_TN, time_TX_TN, covparms, cov_name, ind_Texas_big,
    file = "results/PCE_samp_TN.RData"
  )
}

# Run PCE_validate.R for Figure 6 and PCE_figures.R for Figures 6--7.
