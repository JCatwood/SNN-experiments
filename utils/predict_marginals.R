# Conditional Gaussian marginals, mixed over sampled training responses.
# Residuals are independent across test rows: these are for marginal scores only.
predict_marginals <- function(train_x, test_x, train_draws, covparms,
                              cov_name = "matern15_scaledim", m = 50, seed = 321) {
  kernel <- getFromNamespace(cov_name, "GpGp")
  nn <- RANN::nn2(train_x, test_x, k = min(m, nrow(train_x)))$nn.idx
  mu <- matrix(0, nrow(test_x), ncol(train_draws))
  v <- numeric(nrow(test_x))
  for (i in seq_len(nrow(test_x))) {
    j <- nn[i, ]
    K <- kernel(covparms, rbind(test_x[i, ], train_x[j, , drop = FALSE]))
    w <- solve(K[-1, -1, drop = FALSE], K[-1, 1])
    mu[i, ] <- crossprod(w, train_draws[j, , drop = FALSE])
    v[i] <- max(0, K[1, 1] - sum(w * K[-1, 1]))
  }
  set.seed(seed)
  list(mean = mu, variance = v,
    draws = mu + sqrt(v) * matrix(rnorm(length(mu)), nrow(mu)))
}
