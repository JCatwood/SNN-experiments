library(scoringRules)

score_output <- function(y_test_samp, y_test, comp_time, scene_ID = 1,
                         m = NULL, method = "SNN", parms = "known") {
  y_test_avg <- mean(exp(y_test))
  pred_samp <- colMeans(exp(y_test_samp))
  values <- c(RMSE = sqrt(mean((y_test_avg - pred_samp)^2)),
    CRPS = scoringRules::crps_sample(y_test_avg, pred_samp), time = comp_time)
  result <- data.frame(scenario = scene_ID, replicate = k,
    m = if (is.null(m)) NA_integer_ else m, score = names(values), method = method,
    cov_kernel = parms, value = as.numeric(values))
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  file <- paste0("results/scene_", scene_ID, "_rep_", k, "_", method,
    "_", parms, "_m", if (is.null(m)) "NA" else m, ".csv")
  write.csv(result, file, row.names = FALSE)
  print(result)
}

kriging_score_output <- function(y_samp, y_test, comp_time, scene_ID = 1,
    m = NULL, method = "SNN", parms = "known", train_cov = covmat,
    cross_cov = covmat_train_test, test_cov = covmat_test) {
  W <- solve(train_cov, cross_cov)
  mu <- crossprod(W, y_samp)
  V <- test_cov - crossprod(cross_cov, W)
  V <- (V + t(V)) / 2
  diag(V) <- diag(V) + 1e-10
  set.seed(100000 + k)
  samples <- mu + crossprod(chol(V), matrix(rnorm(length(mu)), nrow(mu)))
  score_output(samples, y_test, comp_time, scene_ID, m, method, parms)
}
