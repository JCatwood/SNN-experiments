library(TruncatedNormal)
library(nntmvn)
library(VeccTMVN)
library(GpGp)
library(R.utils)
library(RANN)
library(ggplot2)
library(tidyr)
library(CensSpBayes)

# input args ---------------
rm(list = ls())
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  scene_ID <- as.integer(args[1])
  n_samp <- as.integer(args[2])
  m <- as.integer(args[3])
  n_burn <- as.integer(args[3])
  thin <- as.integer(args[4])
} else {
  scene_ID <- 1
  n_samp <- 10
  m <- 30
  n_burn <- 2000
  thin <- 5
}

# sim setup ----------------
if (scene_ID == 1) {
  cov_func <- GpGp::matern15_isotropic
  cov_parms <- c(1.0, 0.03, 0.0)
  cov_name <- "matern15_isotropic"
  lb <- -Inf
  ub <- 0
  n_vec <- seq(from = 10, by = 4, length.out = 11)^2
}
n_mtd <- 4
time_df <- matrix(NA, length(n_vec), n_mtd)
max_time <- 1800


# simulation --------------------
n_ind <- 1
for (n_sub in n_vec) {
  cat("n = ", n_sub, "\n")
  set.seed(123)
  if (scene_ID == 1) {
    n_knots <- as.integer(sqrt(n_sub))
    tmp_vec <- seq(from = 0, by = 0.02, length.out = n_knots)
    locs_sub <- as.matrix(expand.grid(tmp_vec, tmp_vec))
    lb_sub <- rep(lb, n_sub)
    ub_sub <- rep(ub, n_sub)
  }
  covmat <- cov_func(cov_parms, locs_sub)
  NN <- nn2(locs_sub, k = m + 1)[[1]]
  # mtd 1 TruncatedNormal
  if (n_ind == 1 || !is.na(time_df[n_ind - 1, 1])) {
    tmp <- R.utils::withTimeout(
      {
        system.time(TruncatedNormal::rtmvnorm(n_samp, rep(0, n_sub),
          sigma = covmat,
          lb = lb_sub, ub = ub_sub
        ))[[3]]
      },
      timeout = max_time,
      onTimeout = "silent",
      elapsed = max_time
    )
    if (is.null(tmp)) {
      time_df[n_ind, 1] <- NA
    } else {
      time_df[n_ind, 1] <- tmp
    }
  } else {
    time_df[n_ind, 1] <- NA
  }

  # mtd 2 VeccTMVN
  if (n_ind == 1 || !is.na(time_df[n_ind - 1, 2])) {
    tmp <- R.utils::withTimeout(
      {
        system.time(VeccTMVN::mvrandn(
          lower = lb_sub, upper = ub_sub, mean = rep(0, n_sub),
          locs = locs_sub, covName = cov_name,
          covParms = cov_parms, m = m, N = n_samp
        ))[[3]]
      },
      timeout = max_time,
      onTimeout = "silent",
      elapsed = max_time
    )
    if (is.null(tmp)) {
      time_df[n_ind, 2] <- NA
    } else {
      time_df[n_ind, 2] <- tmp
    }
  } else {
    time_df[n_ind, 2] <- NA
  }

  # mtd 3 ours
  time_bgn <- Sys.time()
  for (i in 1:n_samp) {
    nntmvn::seq_Vecc_samp_func(
      y = rep(NA, n_sub), cens_lb = lb_sub, cens_ub = ub_sub,
      mask_cens = rep(TRUE, n_sub), NN = NN, covmat = covmat, seed = i
    )
  }
  time_end <- Sys.time()
  time_df[n_ind, 3] <- difftime(time_end, time_bgn, units = "secs")[[1]]

  # mtd 4 CensSpBayes
  time_bgn <- Sys.time()
  inla.mats <- CensSpBayes::create_inla_mats(
    S = locs_sub, S.pred = locs_sub,
    offset = c(0.01, 0.2),
    cutoff = 0.05,
    max.edge = c(0.01, 0.1)
  )
  X.obs <- matrix(1, n_sub, 1)
  X.pred <- matrix(1, n_sub, 1)
  set.seed(123)
  out <- CensSpBayes(
    Y = ub_sub, S = locs_sub, X = X.obs,
    cutoff.Y = ub_sub,
    S.pred = locs_sub, X.pred = X.pred,
    inla.mats = inla.mats,
    rho.init = 0.1, rho.upper = 5,
    iters = n_burn + thin * n_samp, burn = n_burn, thin = thin
  )
  time_end <- Sys.time()
  time_df[n_ind, 4] <- difftime(time_end, time_bgn, units = "secs")[[1]]

  # update misc
  n_ind <- n_ind + 1
}
if (!file.exists("results")) {
  dir.create("results")
}
save(time_df, file = paste0(
  "results/time_for_", n_samp, "_scene", scene_ID, "_m", m,
  ".RData"
))

# plot ------------------------------------------
load(paste0(
  "results/time_for_", n_samp, "_scene", scene_ID, "_m", m,
  ".RData"
))
time_df <- as.data.frame(cbind(n_vec, time_df))
if (sum(is.na(time_df[, 2])) > 0) {
  time_df[which(is.na(time_df[, 2]))[1], 2] <- max_time
}
if (sum(is.na(time_df[, 3])) > 0) {
  time_df[which(is.na(time_df[, 3]))[1], 3] <- max_time
}
if (sum(is.na(time_df[, 4])) > 0) {
  time_df[which(is.na(time_df[, 4]))[1], 4] <- max_time
}
colnames(time_df) <- c("n", "MET", "VMET", "SNN", "CSB")
time_df <- pivot_longer(time_df, c(2:(n_mtd + 1)),
  names_to = "method",
  values_to = "time"
)
ggplot(time_df, aes(x = n, y = time, color = method, linetype = method)) +
  geom_line(size = 2) +
  theme(
    legend.position = c(0.8, 0.8), legend.text = element_text(size = 12),
    legend.title = element_blank(), text = element_text(size = 16)
  ) +
  ylab("seconds")
if (!file.exists("plots")) {
  dir.create("plots")
}
ggsave(paste0(
  "plots/time_for_", n_samp, "_scene", scene_ID, "_m", m,
  ".pdf"
), width = 8, height = 5)
