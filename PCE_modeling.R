library(VeccTMVN)
library(sf)
library(spData)
library(GpGp)
library(CensSpBayes)

rm(list = ls())

run_parm_est <- FALSE
run_VeccTMVN <- FALSE
run_TN <- FALSE
run_CB <- FALSE
run_seq_Vecc <- FALSE
run_CB_all <- FALSE

load("PCE.RData")
summary(data.PCE.censored)
unique(data.PCE.censored$dl_units)
unique(data.PCE.censored$Units)
hist(log(data.PCE.censored$result_va)) # it looks like censored normal indeed!
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
    control = list(trace = 1, maxit = 500), m = 50, NLevel2 = 1e3
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
  ii <- as.integer(st_intersects(pts, states))
  state_names[ii]
}
state_names <- lonlat_to_state(data.frame(locs))
ind_Texas <- which(state_names == "Texas")
# Define a enveloping box for Texas -------------------------------------------
## TX boundary lat: 25.83333 to 36.5 lon: -93.51667 to -106.6333
ind_Texas_box <- which((locs[, 1] < -93.02) & (locs[, 1] > -107.13) &
  (locs[, 2] > 25.33) & (locs[, 2] < 37))
# Sample at locations given by `ind_Texas_big` using VeccTMVN ------------------
if (run_VeccTMVN) {
  load("results/PCE_modeling.RData")
  ind_Texas_big <- ind_Texas
  ind_obs <- which(!is.na(y))
  ind_Texas_big <- union(ind_Texas_big, ind_obs)
  ind_censor_Texas_big <- which(is.na(y[ind_Texas_big]))
  locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
  y_scaled_Texas_big <- y_scaled[ind_Texas_big]
  b_scaled_Texas_big <- b_scaled[ind_Texas_big]
  n_samp <- 1000
  m <- 50
  covparms <- exp(opt_obj$par)
  time_TX_VT <- system.time(
    samp_TX_VT <- ptmvrandn(
      locs_scaled_Texas_big,
      ind_censor_Texas_big,
      y_scaled_Texas_big,
      b_scaled_Texas_big, cov_name, covparms,
      m = m, N = n_samp
    )
  )[[3]]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(samp_TX_VT, time_TX_VT, covparms, cov_name, ind_Texas_big,
    file = "results/PCE_samp_VT.RData"
  )
}
# Seq Vecc method --------------------
if (run_seq_Vecc) {
  library(nntmvn)
  library(RANN)
  library(doParallel)
  load("results/PCE_modeling.RData")
  ind_obs <- which(!is.na(y))
  ind_cens <- which(is.na(y))
  n_samp <- 1000
  n_cens <- sum(is.na(y_scaled))
  covparms <- exp(opt_obj$par)
  m_cens <- 25
  m_obs <- 25
  ## Sample at all locations using nntmvn -----------------
  time_bgn <- Sys.time()
  dim_scale <- diag(c(1, 1, 0.001))
  NN_cens <- RANN::nn2(locs_scaled[ind_cens, ] %*% dim_scale,
    locs_scaled %*% dim_scale,
    k = m_cens + 1
  )[[1]]
  NN_obs <- RANN::nn2(locs_scaled[ind_obs, ] %*% dim_scale,
    locs_scaled %*% dim_scale,
    k = m_obs
  )[[1]]
  NN_cens <- apply(NN_cens, c(1, 2), function(x) {
    ind_cens[x]
  })
  NN_obs <- apply(NN_obs, c(1, 2), function(x) {
    ind_obs[x]
  })
  NN <- cbind(NN_cens, NN_obs)
  all(locs_scaled[NN[ind_cens, 1], ] == locs_scaled[ind_cens, ])
  ncores <- 4
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  samp_seq_Vecc_all <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
    nntmvn::rtmvn_snn(
      y = y_scaled, cens_lb = rep(-Inf, n),
      cens_ub = b_scaled,
      mask_cens = is.na(y_scaled), locs = locs_scaled,
      NN = NN, cov_name = cov_name, cov_parm = covparms
    )
  }
  stopCluster(cl)
  time_end <- Sys.time()
  time_seq_Vecc_all <- difftime(time_end, time_bgn, units = "secs")[[1]]
  samp_seq_Vecc_all <- matrix(unlist(samp_seq_Vecc_all),
    length(samp_seq_Vecc_all[[1]]),
    n_samp,
    byrow = FALSE
  )
  pred_seq_Vecc_all <- rowMeans(samp_seq_Vecc_all)
  sd_seq_Vecc_all <- apply(samp_seq_Vecc_all, 1, sd) / sqrt(n_samp)
  samp_seq_Vecc_all <- samp_seq_Vecc_all[, 1]
  if (!file.exists("results")) {
    dir.create("results")
  }
  save(pred_seq_Vecc_all, sd_seq_Vecc_all, samp_seq_Vecc_all,
    time_seq_Vecc_all, covparms, cov_name,
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
  n_iter_MC <- 25000
  thin <- 5
  ## Sample at locations given by `ind_Texas_big` using CB -----------------
  time_bgn <- Sys.time()
  inla.mats <- create_inla_mats(
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
  y_samp_CB$Y.pred.samp <-
    y_samp_CB$Y.pred.samp[, ncol(y_samp_CB$Y.pred.samp), drop = F]
  save(y_samp_CB, time_TX_CB, ind_Texas_big,
    file = "results/PCE_samp_CB.RData"
  )
}
if (run_CB_all) {
  y_obs_scaled <- y_scaled
  y_obs_scaled[mask_censor] <- b_scaled[mask_censor]
  n_burn <- 20000
  n_iter_MC <- 25000
  thin <- 5
  ## Sample at all locations using CB -----------------
  time_bgn <- Sys.time()
  inla.mats <- create_inla_mats(
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
  y_samp_CB_all$Y.pred.samp <-
    y_samp_CB_all$Y.pred.samp[, ncol(y_samp_CB_all$Y.pred.samp), drop = F]
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
  n_samp <- 1000
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
# Texas plots -------------------------------------------
library(fields)
library(RColorBrewer)
if (!file.exists("plots")) {
  dir.create("plots")
}
lon_grid <- seq(from = -106.6, to = -93.5, length.out = 100)
lat_grid <- seq(from = 25.8, to = 36.6, length.out = 100)
TX_grid <- expand.grid(lon_grid, lat_grid)
colnames(TX_grid) <- c("lon", "lat")
TX_grid_scaled <- TX_grid
for (j in 1:(d - 1)) {
  TX_grid_scaled[, j] <- (TX_grid_scaled[, j] - min(locs[, j])) /
    (max(locs[, j]) - min(locs[, j]))
}
TX_grid_scaled <- cbind(TX_grid_scaled, rep(1, nrow(TX_grid)))
colnames(TX_grid_scaled) <- c("lon", "lat", "time")
## pred VeccTMVN --------------------------------------
load("results/PCE_samp_VT.RData")
locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
cov_mat <- get(cov_name)(covparms, rbind(
  locs_scaled_Texas_big,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_scaled_Texas_big)
n_grid <- nrow(TX_grid)
pred_VMET <- rowMeans(cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], samp_TX_VT))
pred_VMET <- pred_VMET * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_VMET[!(lonlat_to_state(TX_grid) == "Texas") |
  (is.na(lonlat_to_state(TX_grid)))] <- NA
## pred seq Vecc --------------------------------------
load("results/PCE_samp_seq_Vecc_all.RData")
locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
cov_mat <- get(cov_name)(covparms, rbind(
  locs_scaled_Texas_big,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_scaled_Texas_big)
n_grid <- nrow(TX_grid)
pred_seq_Vecc <- as.vector(cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], pred_seq_Vecc_all[ind_Texas_big]))
pred_seq_Vecc <- pred_seq_Vecc * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_seq_Vecc[!(lonlat_to_state(TX_grid) == "Texas") |
  (is.na(lonlat_to_state(TX_grid)))] <- NA
## pred_LOD_GP --------------------------------------
ind_censor_Texas_big <- which(is.na(y[ind_Texas_big]))
locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
b_scaled_Texas_big <- b_scaled[ind_Texas_big]
y_aug <- samp_TX_VT[, 1]
y_aug[ind_censor_Texas_big] <- b_scaled_Texas_big[ind_censor_Texas_big]
cov_mat <- get(cov_name)(covparms, rbind(
  locs_scaled_Texas_big,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_scaled_Texas_big)
n_grid <- nrow(TX_grid)
pred_GP_aug <- as.vector(cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], y_aug))
pred_GP_aug <- pred_GP_aug * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_GP_aug[!(lonlat_to_state(TX_grid) == "Texas") |
  (is.na(lonlat_to_state(TX_grid)))] <- NA
## pred_CB --------------------------------------
load("results/PCE_samp_CB.RData")
ind_censor_Texas_big <- which(is.na(y[ind_Texas_big]))
locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
y_scaled_Texas_big_tmp <- y_scaled[ind_Texas_big]
y_scaled_Texas_big_tmp[is.na(y_scaled_Texas_big_tmp)] <-
  y_samp_CB$Y.pred.posmean
cov_mat <- get(cov_name)(covparms, rbind(
  locs_scaled_Texas_big,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_scaled_Texas_big)
n_grid <- nrow(TX_grid)
pred_GP_CB <- as.vector(cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], y_scaled_Texas_big_tmp))
pred_GP_CB <- pred_GP_CB * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_GP_CB[!(lonlat_to_state(TX_grid) == "Texas") |
  (is.na(lonlat_to_state(TX_grid)))] <- NA
## pred_TN --------------------------------------
load("results/PCE_samp_TN.RData")
ind_censor_Texas_big <- which(is.na(y[ind_Texas_big]))
locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
y_scaled_Texas_big_tmp <- y_scaled[ind_Texas_big]
y_scaled_Texas_big_tmp[is.na(y_scaled_Texas_big_tmp)] <-
  colMeans(samp_TX_TN)
cov_mat <- get(cov_name)(covparms, rbind(
  locs_scaled_Texas_big,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_scaled_Texas_big)
n_grid <- nrow(TX_grid)
pred_TN <- as.vector(cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], y_scaled_Texas_big_tmp))
pred_TN <- pred_TN * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_TN[!(lonlat_to_state(TX_grid) == "Texas") |
  (is.na(lonlat_to_state(TX_grid)))] <- NA
## actual plot ----------------------------------------
zlim <- range(pred_GP_aug, pred_seq_Vecc, pred_GP_CB, pred_TN, pred_VMET,
  na.rm = TRUE
)
plot_TX <- function(pred_vec, mtd) {
  pdf(file = paste0("plots/PCE_TX_", mtd, ".pdf"), width = 5, height = 5)
  image.plot(lon_grid, lat_grid, matrix(pred_vec, 100, 100),
    col = colorRampPalette(brewer.pal(11, "RdBu")[11:1])(30),
    xlab = "longitude", ylab = "latitude", cex.lab = 1.3,
    cex.axis = 1.3, legend.shrink = 0.8, legend.cex = 2.5, legend.width = 2,
    zlim = zlim,
    mgp = c(2, 1, 0)
  )
  points(
    locs[ind_Texas_big[ind_censor_Texas_big], 1],
    locs[ind_Texas_big[ind_censor_Texas_big], 2],
    col = "grey",
    cex = 0.6, pch = 1,
  )
  ind_obs <- which(!is.na(y))
  points(
    x = locs[ind_obs, 1], y = locs[ind_obs, 2], col = "black",
    cex = 0.6, pch = 4,
  )
  dev.off()
}
# fields::US(xlim = c(-106.6, -93.5), ylim = c(25.8, 36.6), add = T)
mtd_vec <- c("GP_aug", "seq_Vecc", "CB", "TN", "VT")
var_names <- c(
  "pred_GP_aug", "pred_seq_Vecc", "pred_GP_CB",
  "pred_TN", "pred_VMET"
)
for (i in 1:5) {
  plot_TX(get0(var_names[i]), mtd_vec[i])
}
# US plot --------------------------------------------
library(autoimage)
load("results/PCE_samp_CB_all.RData")
load("results/PCE_samp_seq_Vecc_all.RData")
zlim <- range(y_samp_CB_all$Y.pred.posmean, samp_seq_Vecc_all) *
  sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)
y_tmp <- y
y_tmp[is.na(y)] <- y_samp_CB_all$Y.pred.posmean * sd(y, na.rm = TRUE) +
  mean(y, na.rm = TRUE)
pdf(file = paste0("plots/PCE_all_CB.pdf"), width = 7, height = 5)
autopoints(locs[, 1], locs[, 2], y_tmp,
  zlim = zlim, col = colorRampPalette(brewer.pal(11, "RdBu")[11:1])(30),
  map = "state", legend = "vertical",
  xlab = "longitude", ylab = "latitude", cex.lab = 1.3,
  cex.axis = 1.3, legend.axis.args = list(at = round(seq(
    from = zlim[1], to = zlim[2], length.out = 7
  ), digits = 2))
)
dev.off()
pdf(file = paste0("plots/PCE_all_seq_Vecc.pdf"), width = 7, height = 5)
y_tmp <- samp_seq_Vecc_all * sd(y, na.rm = TRUE) +
  mean(y, na.rm = TRUE)
autopoints(locs[, 1], locs[, 2], y_tmp,
  zlim = zlim, col = colorRampPalette(brewer.pal(11, "RdBu")[11:1])(30),
  map = "state", legend = "vertical",
  xlab = "longitude", ylab = "latitude", cex.lab = 1.3,
  cex.axis = 1.3, legend.axis.args = list(at = round(seq(
    from = zlim[1], to = zlim[2], length.out = 7
  ), digits = 2))
)
dev.off()
