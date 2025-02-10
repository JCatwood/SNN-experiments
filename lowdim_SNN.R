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
# reorder_mtd
reorder_mtd <- 0
# sim setup ---------------------------------------
tmp_vec <- seq(from = 0.025, to = 0.975, by = 0.05)
locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
rm(tmp_vec)
n <- nrow(locs)
m <- 30
range_parm <- 0.1
nu_parm <- 1.5
levl_censor <- rep(1, n)
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
n_samp <- 50
# generate TMVN samples using multi-threads -------------------
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)
y_samp_SNN <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
  rptmvn(y_obs, rep(-Inf, n), levl_censor, mask_cens,
    m = m, locs = locs,
    covmat = covmat, ordering = reorder_mtd, seed = i
  )
}
stopCluster(cl)
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
cond_covmat_cens[lower.tri(cond_covmat_cens)] <-
  t(cond_covmat_cens)[lower.tri(cond_covmat_cens)]
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
