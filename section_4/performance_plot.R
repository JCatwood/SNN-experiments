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
a <- filter(s, score == !!score, method %in% chosen, cov_kernel == "known")
t <- filter(s, score == "time") %>% select(scenario, m, method, cov_kernel, time = avg_value) %>%
  mutate(time = if_else(!startsWith(method, "SNN"), time + 600, time))
a <- left_join(a, t, by = c("scenario", "m", "method", "cov_kernel"))
a$label <- ifelse(a$method == "VT", "VMET", ifelse(a$method == "MET", "MET",
  ifelse(a$cov_kernel == "unknown", "SNN (fitted)", "SNN")))
a <- a %>% filter(!(label == "VMET" & scenario == "2" & m > 10))
for (scene in sort(unique(a$scenario))) {
  p <- ggplot(filter(a, scenario == scene), aes(time, avg_value, color = label, shape = label)) +
    geom_errorbar(aes(ymin = avg_value - 1.96 * se, ymax = avg_value + 1.96 * se), width = 0) +
    geom_line(aes(group = label)) + geom_point(size = 2.5) +
    geom_text(aes(label = ifelse(is.na(m), "", m)), vjust = -0.8, show.legend = FALSE, size = 3) +
    scale_x_log10() + labs(x = "Time (seconds)", y = score, color = NULL, shape = NULL) +
    theme_bw(base_size = 14) + 
    theme(legend.position = "bottom", legend.text = element_text(size = 14))
  ggsave(paste0("plots/revision/performance_plot_scene_", scene, "_", score, "_order_", order_ID, ".pdf"),
    p, width = 5.5, height = 4.5)
}

method_labels <- c(
  SNN = "No reordering",
  SNN_order_desc = "Descending order",
  SNN_order_maximin = "Maximin order"
)

a <- x %>%
  filter(
    scenario == 1, score == "RMSE", cov_kernel == "known",
    method %in% names(method_labels), !is.na(m)
  ) %>%
  mutate(
    m = factor(m, levels = sort(unique(m))),
    ordering = factor(method_labels[method], levels = unname(method_labels))
  )

if (!nrow(a)) stop("No Scenario 1 RMSE results found for the three SNN orderings.")

p <- ggplot(a, aes(m, value, fill = ordering)) +
  geom_boxplot(width = 0.7, position = position_dodge(width = 0.8), outlier.size = 1.2) +
  scale_fill_manual(values = c("#4E79A7", "#E15759", "#59A14F")) +
  labs(x = "Conditioning size, m", y = "RMSE", fill = "Ordering") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

dir.create("plots/revision", recursive = TRUE, showWarnings = FALSE)
ggsave("plots/revision/order_compare_scene_1_RMSE.pdf", p, width = 7.5, height = 4.5)
print(p)

