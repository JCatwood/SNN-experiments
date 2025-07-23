library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

m_cmp_rslt <- read.table("m_cmp.csv", header = FALSE, sep = ",")
mtd_cmp_rslt <- read.table("mtd_cmp.csv", header = FALSE, sep = ",")
colnames(m_cmp_rslt) <- c("scenario", "m", "score", "method", "cov_kernel", "value")
colnames(mtd_cmp_rslt) <- c("scenario", "score", "method", "cov_kernel", "value")
mtd_cmp_rslt$m <- NA
mtd_cmp_rslt$m[mtd_cmp_rslt$method == "SNN"] <- 30
mtd_cmp_rslt$m[mtd_cmp_rslt$method == "SNN_order_maximin"] <- 30
mtd_cmp_rslt$m[mtd_cmp_rslt$method == "VT"] <- 30
mtd_cmp_rslt <- mtd_cmp_rslt[c(1, 6, 2, 3, 4, 5)]
all_rslt <- rbind(m_cmp_rslt, mtd_cmp_rslt)
all_rslt <- all_rslt[!all_rslt$method == "SNN_order_ascd", ] # remove order_ascd
all_rslt$submethod <- all_rslt$method
ind_SNN <- which(all_rslt$method %in% c("SNN", "SNN_order_desc", "SNN_order_maximin"))
all_rslt$submethod[ind_SNN] <- paste0(
  "SNN",
  "[", all_rslt$m[ind_SNN], "]_",
  all_rslt$cov_kernel[ind_SNN]
)
ind_VT <- which(all_rslt$method == "VT")
all_rslt$submethod[ind_VT] <- paste0(
  "VMET",
  "[", all_rslt$m[ind_VT], "]_",
  all_rslt$cov_kernel[ind_VT]
)
ind_TN <- which(all_rslt$method == "TN")
all_rslt$submethod[ind_TN] <- "MET"
ind_CB <- which(all_rslt$method == "CB")
all_rslt$submethod[ind_CB] <- "CSB"

m_vec <- sort(unique(m_cmp_rslt$m))
m_length <- length(m_vec)

# comparison between SNN and others ---------------------------------------
score_chosen <- "RMSE"
order <- "maximin" # "maximin" # scenario 1 and 2 do not have desc ordering results
if (is.null(order)) {
  mtd_vec <- c("CB", "SNN", "VT", "TN")
} else {
  mtd_vec <- c("CB", paste0("SNN", "_order_", order), "VT", "TN")
}
for (scenario_chosen in sort(unique(m_cmp_rslt$scenario))) {
  subset_score <- all_rslt %>%
    filter(
      score == score_chosen & scenario == scenario_chosen & method %in% mtd_vec
    ) %>%
    group_by(submethod, score) %>%
    summarise(value_avg = mean(value))
  subset_time <- all_rslt %>%
    filter(
      score == "time" & scenario == scenario_chosen & method %in% mtd_vec
    ) %>%
    group_by(submethod, score) %>%
    summarise(value_avg = mean(value))
  subset <- subset_score %>%
    left_join(
      subset_time,
      by = "submethod"
    ) %>%
    rename(score = score.x, value = value_avg.x, time = value_avg.y) %>%
    select(!score.y)

  subset <- subset %>% mutate(submethod = ifelse(
    grepl("SNN.*_known", submethod),
    "SNN", submethod
  ))
  subset <- subset %>% mutate(submethod = ifelse(
    grepl("SNN.*_unknown", submethod),
    "SNN_unknown", submethod
  ))
  subset <- subset %>% mutate(submethod = ifelse(
    grepl("VMET.*_known", submethod),
    "VMET", submethod
  ))
  subset$submethod <- factor(subset$submethod, levels = c(
    "SNN", "SNN_unknown", "VMET", "MET", "CSB"
  ))

  ggplot(data = subset, mapping = aes(x = time, y = value)) +
    geom_point(mapping = aes(colour = submethod, shape = submethod), size = 3) +
    geom_line(mapping = aes(group = submethod, colour = submethod)) +
    scale_y_continuous(
      name = score_chosen
    ) +
    scale_x_continuous(
      name = "time (seconds)", trans = "log2"
    ) +
    scale_color_manual(labels = c(
      expression("SNN"),
      expression(SNN^"*"),
      expression("VMET"),
      expression("MET"),
      expression(CSB^"*")
    ), values = hue_pal()(5)) +
    scale_shape_manual(values = 1:5, labels = c(
      expression("SNN"),
      expression(SNN^"*"),
      expression("VMET"),
      expression("MET"),
      expression(CSB^"*")
    )) +
    # ggtitle(paste("Scenario", prob_ind)) +
    theme(
      text = element_text(size = 16),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.13, 0.8)
    )
  if (is.null(order)) {
    order_ID <- 0
  } else if (order == "desc") {
    order_ID <- 1
  } else if (order == "maximin") {
    order_ID <- 2
  } else {
    stop("Wrong order name\n")
  }
  ggsave(
    paste0(
      "plots/performance_plot_scene_", scenario_chosen, "_", score_chosen,
      "_order_", order_ID, ".pdf"
    ),
    width = 5.5,
    height = 5
  )
}

# comparison of different orderings using SNN --------------------------------
score_chosen <- "CRPS"
mtd_vec <- c("SNN", "SNN_order_desc", "SNN_order_maximin")
scenario_chosen <- 3

subset_mask <- all_rslt$score == score_chosen & all_rslt$scenario == scenario_chosen &
  all_rslt$method %in% mtd_vec
subset_rslt <- all_rslt[subset_mask, , drop = FALSE]
subset_rslt$method <- factor(subset_rslt$method,
  levels = c("SNN", "SNN_order_desc", "SNN_order_maximin")
)
subset_rslt$m <- factor(subset_rslt$m, levels = sort(unique(subset_rslt$m)))
ggplot(data = subset_rslt, mapping = aes(x = m, y = value)) +
  geom_boxplot(mapping = aes(fill = method), alpha = 0.5) +
  scale_y_continuous(name = score_chosen) +
  scale_x_discrete(labels = paste0("m = ", unique(subset_rslt$m))) +
  scale_fill_discrete(labels = c("default", "var-desc", "maximin")) +
  theme(
    legend.title = element_blank(),
    text = element_text(size = 16),
    legend.position = c(0.8, 0.8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
ggsave(
  paste0("plots/order_compare_scene_", scenario_chosen, "_", score_chosen, ".pdf"),
  width = 5,
  height = 5
)
