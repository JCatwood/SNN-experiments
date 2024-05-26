library(fields)
library(CensSpBayes)

# sim setup ---------------------------------------
rm(list = ls())
set.seed(123)
tmp_vec <- seq(from = 0.025, to = 0.975, by = 0.05)
locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
tmp_vec <- seq(from = 0, to = 1.0, by = 0.02)
rm(tmp_vec)
n <- nrow(locs)
range_parm <- 0.1
nu_parm <- 1.5
covmat <- fields::Matern(as.matrix(dist(locs)), range = range_parm, nu = nu_parm)
L <- t(chol(covmat))
y <- as.vector(L %*% rnorm(n))
levl_censor <- rep(1.0, n)
mask_cens <- y < levl_censor
ind_cens <- which(mask_cens)
n_cens <- length(ind_cens)
burn <- 2000
iter <- 2250
thin <- 5
n_samp <- floor((iter - burn) / thin)

y_obs <- y
y_obs[mask_cens] <- levl_censor[mask_cens]

# CensSpBayes -----------------------------------
inla.mats <- create_inla_mats(
  S = locs, S.pred = locs[mask_cens, ],
  offset = c(0.01, 0.2),
  cutoff = 0.05,
  max.edge = c(0.01, 0.1)
)

X.obs <- matrix(1, nrow(locs), 1)
X.pred <- matrix(1, length(ind_cens), 1)
set.seed(123)
out <- CensSpBayes(
  Y = y_obs, S = locs, X = X.obs,
  cutoff.Y = levl_censor,
  S.pred = locs[mask_cens, ], X.pred = X.pred,
  inla.mats = inla.mats,
  rho.init = 0.1, rho.upper = 5,
  iters = iter, burn = burn, thin = thin, ret_samp = TRUE
)
samp_CSB <- out$Y.pred.samp
if (!file.exists("plots")) {
  dir.create("plots")
}
# sample TMVN with TruncatedNormal ---------------------------------------
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
# qq plot -------------------------
pdf(
  file = paste0("plots/lowdim_scene1_CSB.pdf"),
  width = 5, height = 5
)
par(mar = c(4, 4, 1, 1))
qqplot(as.vector(y_cens_samp_MET), as.vector(samp_CSB),
  xlab = "samples from MET",
  ylab = "samples from CSB",
  xlim = range(y_cens_samp_MET, samp_CSB),
  ylim = range(y_cens_samp_MET, samp_CSB),
  cex.lab = 1.3, cex.axis = 1.3
)
abline(a = 0, b = 1, col = "red", lwd = 2)
dev.off()
# compute RMSE ------------------------------------
cat(
  "RMSE for CSB is ",
  sqrt(mean((out$Y.pred.posmean - y[mask_cens])^2)), "\n"
)
cat(
  "RMSE for MET is ",
  sqrt(mean((rowMeans(y_cens_samp_MET) - y[mask_cens])^2)), "\n"
)

# compute nll ------------------------------------
cat("nll for CSB is ", -mean(
  dnorm(y[mask_cens],
    mean = out$Y.pred.posmean,
    sd = sqrt(out$Y.pred.posvar) / sqrt(n_samp), log = T
  )
), "\n")
cat("nll for MET is ", -mean(
  dnorm(y[mask_cens],
    mean = rowMeans(y_cens_samp_MET),
    sd = apply(y_cens_samp_MET, 1, sd) / sqrt(n_samp), log = T
  )
), "\n")
