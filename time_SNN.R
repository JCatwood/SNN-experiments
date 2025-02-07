library(TruncatedNormal)
library(nntmvn)
library(VeccTMVN)
library(GpGp)
library(R.utils)
library(RANN)
library(ggplot2)
library(tidyr)

# input args ---------------
rm(list = ls())
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  scene_ID <- as.integer(args[1])
  n_samp <- as.integer(args[2])
  m <- as.integer(args[3])
} else {
  scene_ID <- 1
  n_samp <- 10
  m <- 30
}

# sim setup ----------------
if (scene_ID == 1) {
  cov_func <- GpGp::matern15_isotropic
  cov_parms <- c(1.0, 0.03, 0.0)
  cov_name <- "matern15_isotropic"
  lb <- -Inf
  ub <- 0
  n_vec <- seq(from = 50, by = 50, length.out = 6)^2
}
n_mtd <- 1
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
  NN <- nn2(locs_sub, k = m + 1)[[1]]
  # mtd 1 ours
  time_bgn <- Sys.time()
  for (i in 1:n_samp) {
    nntmvn::rptmvn(
      y = rep(NA, n_sub), cens_lb = lb_sub, cens_ub = ub_sub,
      mask_cens = rep(TRUE, n_sub), NN = NN, cov_name = cov_name, locs = locs_sub,
      cov_parm = cov_parms, seed = i
    )
  }
  time_end <- Sys.time()
  time_df[n_ind, 1] <- difftime(time_end, time_bgn, units = "secs")[[1]]

  # update misc
  n_ind <- n_ind + 1
}
if (!file.exists("results")) {
  dir.create("results")
}
save(time_df, file = paste0(
  "results/timeSNN_for_", n_samp, "_scene", scene_ID, "_m", m,
  ".RData"
))

# plot ------------------------------------------
load(paste0(
  "results/timeSNN_for_", n_samp, "_scene", scene_ID, "_m", m,
  ".RData"
))
time_df <- as.data.frame(cbind(n_vec, time_df))
if (sum(is.na(time_df[, 2])) > 0) {
  time_df[which(is.na(time_df[, 2]))[1], 2] <- max_time
}
colnames(time_df) <- c("n", "SNN")
time_df <- pivot_longer(time_df, c(2:(n_mtd + 1)),
  names_to = "method",
  values_to = "time"
)
ggplot(time_df, aes(x = n, y = time)) +
  geom_line(linewidth = 2, color = "#00BFC4", linetype = "longdash") +
  theme(
    legend.position.inside = c(0.8, 0.8), legend.text = element_text(size = 12),
    legend.title = element_blank(), text = element_text(size = 16)
  ) +
  ylab("seconds")
if (!file.exists("plots")) {
  dir.create("plots")
}
ggsave(paste0(
  "plots/timeSNN_for_", n_samp, "_scene", scene_ID, "_m", m,
  ".pdf"
), width = 8, height = 5)
