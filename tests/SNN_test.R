library(fields)
library(TruncatedNormal)
library(mvtnorm)
library(RANN)
library(doParallel)
library(GpGp)
library(nntmvn)
library(scoringRules)

# input args ---------------------------------------
rm(list = ls())
set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  reorder_mtd <- as.numeric(args[1])
} else {
  # reorder_mtd
  #  0 for no reordering
  #  1 for maximin
  #  2 for Vecchia as in VeccTMVN
  reorder_mtd <- 0
}
# sim setup ---------------------------------------
tmp_vec <- seq(from = 0.025, to = 0.975, by = 0.05)
locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
rm(tmp_vec)
n <- nrow(locs)
m <- 30
range_parm <- 0.1
nu_parm <- 1.5
levl_censor <- rep(1.0, n)
covmat <- fields::Matern(as.matrix(dist(locs)),
  range = range_parm,
  nu = nu_parm
)
L <- t(chol(covmat))
y <- as.vector(L %*% rnorm(n))
mask_cens <- y < levl_censor
n_cens <- sum(mask_cens)
y_obs <- y
y_obs[mask_cens] <- NA # levl_censor[mask_cens]
NN <- nn2(locs, k = m + 1)[[1]]
n_samp <- 50
all(NN[, 1] == 1:n)
order_new <- NULL
if (reorder_mtd == 1) { # maximin
  order_new <- GpGp::order_maxmin(locs)
  locs_order <- locs[order_new, , drop = FALSE]
  covmat <- fields::Matern(as.matrix(dist(locs_order)),
    range = range_parm,
    nu = nu_parm
  )
  y_obs <- y_obs[order_new]
  mask_cens <- mask_cens[order_new]
  NN <- nn2(locs_order, k = m + 1)[[1]]
  levl_censor <- levl_censor[order_new]
} else if (reorder_mtd == 2) { # Vecc
  order_new <- VeccTMVN::Vecc_reorder(rep(-Inf, n), y_obs, m, covMat = covmat)$order
  locs_order <- locs[order_new, , drop = FALSE]
  covmat <- fields::Matern(as.matrix(dist(locs_order)),
    range = range_parm,
    nu = nu_parm
  )
  y_obs <- y_obs[order_new]
  mask_cens <- mask_cens[order_new]
  NN <- nn2(locs_order, k = m + 1)[[1]]
  levl_censor <- levl_censor[order_new]
}
# generate TMVN samples using multi-threads -------------------
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)
y_samp_SNN <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
  rtmvn_snn(y_obs, rep(-Inf, n), levl_censor, mask_cens,
    m = m,
    covmat = covmat, seed = i
  )
}
stopCluster(cl)
if (!is.null(order_new)) {
  order_new_rev <- 1:n
  order_new_rev[order_new] <- 1:n
  y_samp_SNN <- lapply(y_samp_SNN, function(x) {
    x[order_new_rev]
  })
  mask_cens <- mask_cens[order_new_rev]
  NN <- nn2(locs, k = m + 1)[[1]]
  levl_censor <- levl_censor[order_new_rev]
}
y_samp_SNN <- matrix(unlist(y_samp_SNN), n, n_samp, byrow = FALSE)
y_cens_samp_SNN <- y_samp_SNN[mask_cens, ]

# sample TMVN with TruncatedNormal ---------------------------------------
covmat <- fields::Matern(as.matrix(dist(locs)),
  range = range_parm,
  nu = nu_parm
)
cond_mean_cens <- as.vector(covmat[mask_cens, !mask_cens] %*%
  solve(covmat[!mask_cens, !mask_cens]) %*%
  y[!mask_cens])
cond_covmat_cens <- covmat[mask_cens, mask_cens] - covmat[mask_cens, !mask_cens] %*%
  solve(covmat[!mask_cens, !mask_cens]) %*% covmat[!mask_cens, mask_cens]
y_cens_samp_MET <- t(TruncatedNormal::rtmvnorm(n_samp,
  mu = cond_mean_cens, sigma = cond_covmat_cens,
  lb = rep(-Inf, n_cens),
  ub = levl_censor[mask_cens]
))
plot(y[mask_cens], rowMeans(y_cens_samp_MET),
  xlab = "true y at censored locs",
  ylab = "predicted y by MET",
  xlim = range(y[mask_cens], rowMeans(y_cens_samp_MET)),
  ylim = range(y[mask_cens], rowMeans(y_cens_samp_MET))
)
abline(a = 0, b = 1, col = "red")
plot(y[mask_cens], rowMeans(y_cens_samp_SNN),
  xlab = "true y at censored locs",
  ylab = "predicted y by SNN",
  xlim = range(y[mask_cens], rowMeans(y_cens_samp_SNN)),
  ylim = range(y[mask_cens], rowMeans(y_cens_samp_SNN))
)
abline(a = 0, b = 1, col = "red")
if (!file.exists("plots")) {
  dir.create("plots")
}
pdf(
  file = paste0("plots/lowdim_scene1_seq_Vecc_order", reorder_mtd, ".pdf"),
  width = 5, height = 5
)
par(mar = c(4, 4, 1, 1))
qqplot(as.vector(y_cens_samp_MET), as.vector(y_cens_samp_SNN),
  xlab = "samples from MET",
  ylab = "samples from SNN",
  xlim = range(y_cens_samp_MET, y_cens_samp_SNN),
  ylim = range(y_cens_samp_MET, y_cens_samp_SNN),
  cex.lab = 1.3, cex.axis = 1.3
)
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()

# compute log-score for the samples -----------------------------
mean(mvtnorm::dmvnorm(t(y_cens_samp_SNN),
  mean = cond_mean_cens,
  sigma = cond_covmat_cens, log = TRUE
))
sd(mvtnorm::dmvnorm(t(y_cens_samp_SNN),
  mean = cond_mean_cens,
  sigma = cond_covmat_cens, log = TRUE
)) / sqrt(n_samp)
mean(mvtnorm::dmvnorm(t(y_cens_samp_MET),
  mean = cond_mean_cens,
  sigma = cond_covmat_cens, log = TRUE
))
sd(mvtnorm::dmvnorm(t(y_cens_samp_MET),
  mean = cond_mean_cens,
  sigma = cond_covmat_cens, log = TRUE
)) / sqrt(n_samp)
cat("Log-score may not make sense as we are not estimating a model\n")

# compute RMSE ------------------------------------
cat(
  "RMSE for SNN is ",
  sqrt(mean((rowMeans(y_cens_samp_SNN) - y[mask_cens])^2)), "\n"
)
cat(
  "RMSE for MET is ",
  sqrt(mean((rowMeans(y_cens_samp_MET) - y[mask_cens])^2)), "\n"
)

# compute nll ------------------------------------
cat("nll for SNN is ", -mean(
  dnorm(y[mask_cens],
    mean = rowMeans(y_cens_samp_SNN),
    sd = apply(y_cens_samp_SNN, 1, sd) / sqrt(n_samp), log = T
  )
), "\n")
cat("nll for MET is ", -mean(
  dnorm(y[mask_cens],
    mean = rowMeans(y_cens_samp_MET),
    sd = apply(y_cens_samp_MET, 1, sd) / sqrt(n_samp), log = T
  )
), "\n")

# compute CRPS ------------------------------------
cat("CRPS for SNN is ", mean(
  scoringRules::crps_sample(y = y[mask_cens], dat = y_cens_samp_SNN)
), "\n")
cat("CRPS for MET is ", mean(
  scoringRules::crps_sample(y = y[mask_cens], dat = y_cens_samp_MET)
), "\n")
