library(ggplot2)
dir.create("plots", showWarnings = FALSE)
# Figure 6: identical held-out splits; points show individual split results.
files <- list.files("results", "^PCE_validation_(TX|US)_[0-9]+\\.csv$", full.names = TRUE)
if (length(files)) {
  scores <- do.call(rbind, lapply(files, read.csv))
  p <- ggplot(subset(scores, metric == "Brier"), aes(method, value, color = method)) +
    geom_boxplot() +
    facet_wrap(~ region, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = "Brier score") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  ggsave("plots/PCE_validation.pdf", p, width = 5, height = 3)
  write.csv(scores, "results/PCE_validation_summary.csv", row.names = FALSE)
}
# Figure 7: posterior mean versus posterior mean, SD versus SD, at censored sites.
load("PCE.RData")
dat <- data.PCE.censored
cens <- as.logical(dat$left_censored)
y <- log(dat$result_va + 1e-8)
center <- mean(y, na.rm = TRUE)
scale_y <- sd(y, na.rm = TRUE)
paths <- c(SNN = "results/PCE_samp_seq_Vecc_all.RData", CSB = "results/PCE_samp_CB_all.RData")
summary <- list()
for (method in names(paths)[file.exists(paths)]) {
  e <- new.env()
  load(paths[method], envir = e)
  X <- if (method == "SNN") e$samp_seq_Vecc_all else e$y_samp_CB_all$Y.pred.samp
  if (!is.matrix(X) || ncol(X) < 2) stop("Regenerate national draws; a single draw cannot estimate SD.")
  if (method == "SNN") X <- X[cens, , drop = FALSE]
  base <- data.frame(lon = dat$lon[cens], lat = dat$lat[cens], method = method)
  summary[[method]] <- rbind(transform(base, statistic = "Posterior mean",
    value = rowMeans(X) * scale_y + center), transform(base, statistic = "Posterior SD",
    value = apply(X, 1, sd) * scale_y))
}
if (length(summary)) {
  maps <- do.call(rbind, summary)
  # Separate files permit one common scale per statistic across all methods.
  for (stat in unique(maps$statistic)) {
    p <- ggplot(subset(maps, statistic == stat), aes(lon, lat, color = value)) +
      geom_point(size = .35, alpha = .8) + 
      # facet_wrap(~ method, nrow = 1) +
      coord_quickmap() + scale_color_viridis_c(name = paste(stat, "(log concentration)")) +
      labs(x = "Longitude", y = "Latitude") + theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
    ggsave(paste0("plots/PCE_national_", if (stat == "Posterior mean") "mean" else "sd", ".pdf"),
      p, width = 5, height = 4)
  }
  write.csv(maps, "results/PCE_national_summary.csv", row.names = FALSE)
}
