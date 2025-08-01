library(scoringRules)

score_output <- function(y_test_samp, y_test, comp_time, scene_ID = c(1, 2, 3),
                         m = NULL, method = c("SNN", "VMET", "MET", "CSB"),
                         parms = c("known", "unknown")) {
  scene_ID <- scene_ID[1]
  method <- method[1]
  parms <- parms[1]

  y_pred <- rowMeans(y_test_samp)

  if (is.null(m)) {
    cat(
      "> ", scene_ID, ", RMSE,", paste(method, parms, sep = ", "), ",",
      sqrt(mean((y_test - y_pred)^2)), "\n"
    )
    cat(
      "> ", scene_ID, ", CRPS,", paste(method, parms, sep = ", "), ",",
      mean(scoringRules::crps_sample(
        y = y_test,
        dat = y_test_samp
      )), "\n"
    )
    cat(
      "> ", scene_ID, ", time,", paste(method, parms, sep = ", "), ",",
      comp_time, "\n"
    )
  } else {
    cat(
      "> ", scene_ID, ",", m, ", RMSE,", paste(method, parms, sep = ", "), ",",
      sqrt(mean((y_test - y_pred)^2)), "\n"
    )
    cat(
      "> ", scene_ID, ",", m, ", CRPS,", paste(method, parms, sep = ", "), ",",
      mean(scoringRules::crps_sample(
        y = y_test,
        dat = y_test_samp
      )), "\n"
    )
    cat(
      "> ", scene_ID, ",", m, ", time,", paste(method, parms, sep = ", "), ",",
      comp_time, "\n"
    )
  }
}
