#####################################################################
#
# Produces the MCMC diagnosis plots for CensSpBayes
#
#####################################################################

library(ggplot2)
library(fields)
library(CensSpBayes)
library(scoringRules)

rm(list = ls())
set.seed(123)

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
for (i in 1:6) {
  mydf <- data.frame(
    x = 1:ncol(samp_CSB),
    y = samp_CSB[i, ]
  )
  ggplot(data = mydf, mapping = aes(x = x, y = y)) +
    geom_point(size = 1) +
    scale_x_continuous(name = "MCMC sample index after burn-in and thinning") +
    scale_y_continuous(name = paste0(i, "-th censored response")) +
    theme(
      text = element_text(size = 14),
      plot.title = element_blank()
    )
  if (!file.exists("plots")) {
    dir.create("plots")
  }
  ggsave(
    paste0("plots/CSB_diagnosis_", i, ".pdf"),
    width = 5,
    height = 5
  )
}
