# command-line arguments ---------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage:",
  "  Rscript performance_plot.R [score] [order]",
  "  Rscript performance_plot.R --score=RMSE --order=maximin",
  "",
  "score: RMSE or CRPS (default: RMSE)",
  "order: none, maximin, or desc (default: none)",
  sep = "\n"
)

if (any(args %in% c("-h", "--help"))) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

score_chosen <- "RMSE"
order_arg <- "none"
if (length(args) > 0) {
  if (any(startsWith(args, "--"))) {
    if (!all(grepl("^--(score|order)=.+$", args))) {
      stop(usage, call. = FALSE)
    }
    for (arg in args) {
      if (startsWith(arg, "--score=")) {
        score_chosen <- sub("^--score=", "", arg)
      } else if (startsWith(arg, "--order=")) {
        order_arg <- sub("^--order=", "", arg)
      }
    }
  } else {
    if (length(args) > 2) {
      stop(usage, call. = FALSE)
    }
    score_chosen <- args[1]
    if (length(args) == 2) {
      order_arg <- args[2]
    }
  }
}

score_chosen <- toupper(score_chosen)
order_arg <- tolower(order_arg)
if (!score_chosen %in% c("RMSE", "CRPS")) {
  stop("score must be RMSE or CRPS.\n", usage, call. = FALSE)
}
if (!order_arg %in% c("none", "null", "na", "maximin", "desc")) {
  stop("order must be none, maximin, or desc.\n", usage, call. = FALSE)
}
order_chosen <- if (order_arg %in% c("none", "null", "na")) NULL else order_arg

library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

if (!dir.exists("plots")) {
  dir.create("plots")
}

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
# Scenarios 1 and 2 do not have variance-descending ordering results.
if (is.null(order_chosen)) {
  mtd_vec <- c("CB", "SNN", "VT", "TN")
} else {
  mtd_vec <- c("CB", paste0("SNN", "_order_", order_chosen), "VT", "TN")
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
  if (is.null(order_chosen)) {
    order_ID <- 0
  } else if (order_chosen == "desc") {
    order_ID <- 1
  } else if (order_chosen == "maximin") {
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
