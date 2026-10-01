library(ggplot2)
library(dplyr)
library(tidyr)
args <- commandArgs(TRUE)
score <- if (length(args)) sub("^--score=", "", args[1]) else "RMSE"
ordering <- if (length(args) > 1) sub("^--order=", "", args[2]) else "maximin"
order_ID <- switch(ordering, none = 0, maximin = 2, desc = 1)
stopifnot(score %in% c("RMSE", "CRPS"), !is.null(order_ID))
files <- list.files("results", "\\.csv$", full.names = TRUE)
if (!length(files)) stop("Run the revised experiments first; legacy CSVs have incompatible scores/designs.")
x <- bind_rows(lapply(files, read.csv))
# Replicates are explicit; repeated jobs overwrite their own result file.
s <- x %>% group_by(scenario, m, method, cov_kernel, score) %>%
  summarise(avg_value = mean(value), se = sd(value) / sqrt(n()), replicates = n(), .groups = "drop")
# write.csv(filter(s, method == "CB"), "results/CSB_summary.csv", row.names = FALSE)
dir.create("plots/revision", recursive = TRUE, showWarnings = FALSE)
chosen <- c(if (ordering == "none") "SNN" else paste0("SNN_order_", ordering), "VT", "MET")
a <- filter(s, score == !!score, method %in% chosen)
t <- filter(s, score == "time") %>% select(scenario, m, method, cov_kernel, time = avg_value)
a <- left_join(a, t, by = c("scenario", "m", "method", "cov_kernel"))
a$label <- ifelse(a$method == "VT", "VMET", ifelse(a$method == "MET", "MET",
  ifelse(a$cov_kernel == "unknown", "SNN (fitted)", "SNN")))
for (scene in sort(unique(a$scenario))) {
  p <- ggplot(filter(a, scenario == scene), aes(time, avg_value, color = label, shape = label)) +
    geom_errorbar(aes(ymin = avg_value - 1.96 * se, ymax = avg_value + 1.96 * se), width = 0) +
    geom_line(aes(group = label)) + geom_point(size = 2.5) +
    geom_text(aes(label = ifelse(is.na(m), "", m)), vjust = -0.8, show.legend = FALSE, size = 3) +
    scale_x_log10() + labs(x = "Time (seconds)", y = score, color = NULL, shape = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(paste0("plots/revision/performance_plot_scene_", scene, "_", score, "_order_", order_ID, ".pdf"),
    p, width = 5.5, height = 4.5)
}
# write.csv(s, "results/performance_summary.csv", row.names = FALSE)
# Ordering comparison is still restricted to Scenario 3.
# p <- ggplot(filter(x, scenario == 3, score == !!score, grepl("^SNN", method)),
#   aes(factor(m), value, fill = method)) + geom_boxplot() + facet_wrap(~ cov_kernel) +
#   labs(x = "Neighborhood size m", y = score, fill = NULL) + theme_bw(base_size = 11) +
#   theme(legend.position = "bottom")
# ggsave(paste0("plots/revision/order_compare_scene_3_", score, ".pdf"), p, width = 8, height = 4)
