m_cmp_rslt <- read.table("m_cmp.csv", header = FALSE, sep = ",")
mtd_cmp_rslt <- read.table("mtd_cmp.csv", header = FALSE, sep = ",")
colnames(m_cmp_rslt) <- c("scenario", "m", "score", "method", "cov_kernel", "value")
colnames(mtd_cmp_rslt) <- c("scenario", "score", "method", "cov_kernel", "value")
mtd_cmp_rslt$m <- NA
mtd_cmp_rslt$m[mtd_cmp_rslt$method == "SNN"] <- 30
mtd_cmp_rslt <- mtd_cmp_rslt[c(1, 6, 2, 3, 4, 5)]
all_rslt <- rbind(m_cmp_rslt, mtd_cmp_rslt)

score <- "CRPS"
m_vec <- sort(unique(m_cmp_rslt$m))
m_length <- length(m_vec)
for (scenario in sort(unique(m_cmp_rslt$scenario))) {
  score_vec <- m_length + 4
  score_sd_vec <- m_length + 4
  ind <- 1
  for (cov_kernel in c("known", "unknown")) {
    if (cov_kernel == "known") {
      for (method in c("VT", "TN", "SNN")) {
        if (method == "SNN") {
          for (m in m_vec) {
            mask <- all_rslt$scenario == scenario & all_rslt$m == m &
              all_rslt$score == score & all_rslt$cov_kernel == cov_kernel
            mask[is.na(mask)] <- FALSE
            score_vec[ind] <- mean(all_rslt$value[mask])
            score_sd_vec[ind] <- sd(all_rslt$value[mask]) / sqrt(sum(mask))
            ind <- ind + 1
          }
        } else {
          mask <- all_rslt$scenario == scenario & all_rslt$method == method &
            all_rslt$score == score & all_rslt$cov_kernel == cov_kernel
          score_vec[ind] <- mean(all_rslt$value[mask])
          score_sd_vec[ind] <- sd(all_rslt$value[mask]) / sqrt(sum(mask))
          ind <- ind + 1
        }
      }
    } else {
      for (method in c("CB", "SNN")) {
        mask <- all_rslt$scenario == scenario & all_rslt$method == method &
          all_rslt$score == score & all_rslt$cov_kernel == cov_kernel
        score_vec[ind] <- mean(all_rslt$value[mask])
        score_sd_vec[ind] <- sd(all_rslt$value[mask]) / sqrt(sum(mask))
        ind <- ind + 1
      }
    }
  }
  cat(sprintf("\\multirow{2}{*}{Scene %d}", scenario),
    sapply(score_vec, FUN = function(x) {
      sprintf("%.2f", x)
    }),
    sep = " & "
  )
  cat("\\\\\n")
  cat(" ", sapply(score_sd_vec, FUN = function(x) {
    sprintf("(%.3f)", x)
  }), sep = " & ")
  cat("\\\\\n")
  cat("\\hline\n")
}
