library(scoringRules)

score_output <- function(y_cens_samp, y_cens, comp_time, scene_ID = 1,
                         m = NULL, method = "SNN", parms = "known") {
  y_cens_pred <- rowMeans(y_cens_samp)
  values <- c(RMSE = sqrt(mean((y_cens - y_cens_pred)^2)),
              CRPS = mean(scoringRules::crps_sample(y_cens, y_cens_samp)), time = comp_time)
  result <- data.frame(scenario = scene_ID, replicate = k,
                       m = if (is.null(m)) NA_integer_ else m, score = names(values), method = method,
                       cov_kernel = parms, value = as.numeric(values))
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  file <- paste0("results/scene_", scene_ID, "_rep_", k, "_", method,
                 "_", parms, "_m", if (is.null(m)) "NA" else m, ".csv")
  write.csv(result, file, row.names = FALSE)
  print(result)
}
