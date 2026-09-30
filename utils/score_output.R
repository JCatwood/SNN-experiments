library(scoringRules)

score_output <- function(y_test_samp, y_test, comp_time, scene_ID = 1,
                         m = NULL, method = "SNN", parms = "known",
                         pred_mean = rowMeans(y_test_samp)) {
  values <- c(RMSE = sqrt(mean((y_test - pred_mean)^2)),
              CRPS = mean(scoringRules::crps_sample(y_test, y_test_samp)), time = comp_time)
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
                                 cross_cov = covmat_train_test, test_variance = diag(covmat_test)) {
  W <- solve(train_cov, cross_cov)
  mu <- crossprod(W, y_samp)
  v <- pmax(0, test_variance - colSums(cross_cov * W))
  # Independent row residuals suffice for marginal CRPS; not joint test draws.
  set.seed(100000 + k)
  samples <- mu + sqrt(v) * matrix(rnorm(length(mu)), nrow(mu))
  score_output(samples, y_test, comp_time, scene_ID, m, method, parms, rowMeans(mu))
}