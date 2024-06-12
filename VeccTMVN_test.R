library(fields)
library(TruncatedNormal)
library(GpGp)
library(VeccTMVN)
library(scoringRules)

# sim setup ---------------------------------------
rm(list = ls())
set.seed(123)
tmp_vec <- seq(from = 0.025, to = 0.975, by = 0.05)
locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
rm(tmp_vec)
n <- nrow(locs)
m <- 50
range_parm <- 0.1
nu_parm <- 1.5
covmat <- fields::Matern(as.matrix(dist(locs)), range = range_parm, nu = nu_parm)
L <- t(chol(covmat))
y <- as.vector(L %*% rnorm(n))
levl_censor <- rep(1.0, n)
mask_cens <- y < levl_censor
ind_cens <- which(mask_cens)
n_cens <- sum(mask_cens)
y_obs <- y
y_obs[mask_cens] <- levl_censor[mask_cens]
n_samp <- 50
# sample with VeccTMVN ----------------------------
samp_VeccTMVN_all <- ptmvrandn(locs, ind_cens, y_obs, levl_censor,
  "matern15_isotropic", c(1.0, range_parm, 0.0),
  m = m, N = n_samp
)
samp_VeccTMVN <- samp_VeccTMVN_all[mask_cens, ]
# sample TMVN with TruncatedNormal ---------------------------------------
cond_mean_cens <- as.vector(covmat[mask_cens, !mask_cens] %*%
  solve(covmat[!mask_cens, !mask_cens]) %*%
  y[!mask_cens])
cond_covmat_cens <- covmat[mask_cens, mask_cens] - covmat[mask_cens, !mask_cens] %*%
  solve(covmat[!mask_cens, !mask_cens]) %*% covmat[!mask_cens, mask_cens]
y_cens_samp_MET <- t(TruncatedNormal::rtmvnorm(50,
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
plot(y[mask_cens], rowMeans(samp_VeccTMVN),
  xlab = "true y at censored locs",
  ylab = "predicted y by Seq Vecc",
  xlim = range(y[mask_cens], rowMeans(samp_VeccTMVN)),
  ylim = range(y[mask_cens], rowMeans(samp_VeccTMVN))
)
abline(a = 0, b = 1, col = "red")
if (!file.exists("plots")) {
  dir.create("plots")
}
pdf(
  file = paste0("plots/lowdim_scene1_VeccTMVN.pdf"),
  width = 5, height = 5
)
par(mar = c(4, 4, 1, 1))
qqplot(as.vector(y_cens_samp_MET), as.vector(samp_VeccTMVN),
  xlab = "samples from MET",
  ylab = "samples from VMET",
  xlim = range(y_cens_samp_MET, samp_VeccTMVN),
  ylim = range(y_cens_samp_MET, samp_VeccTMVN),
  cex.lab = 1.3, cex.axis = 1.3
)
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()
# compute RMSE ------------------------------------
cat(
  "RMSE for VMET is ",
  sqrt(mean((rowMeans(samp_VeccTMVN) - y[mask_cens])^2)), "\n"
)
cat(
  "RMSE for MET is ",
  sqrt(mean((rowMeans(y_cens_samp_MET) - y[mask_cens])^2)), "\n"
)

# compute nll ------------------------------------
cat("nll for VMET is ", -mean(
  dnorm(y[mask_cens],
    mean = rowMeans(samp_VeccTMVN),
    sd = apply(samp_VeccTMVN, 1, sd) / sqrt(n_samp), log = T
  )
), "\n")
cat("nll for MET is ", -mean(
  dnorm(y[mask_cens],
    mean = rowMeans(y_cens_samp_MET),
    sd = apply(y_cens_samp_MET, 1, sd) / sqrt(n_samp), log = T
  )
), "\n")

# compute CRPS ------------------------------------
cat("CRPS for VMET is ", mean(
  scoringRules::crps_sample(y = y[mask_cens], dat = samp_VeccTMVN)
), "\n")
cat("CRPS for MET is ", mean(
  scoringRules::crps_sample(y = y[mask_cens], dat = y_cens_samp_MET)
), "\n")

