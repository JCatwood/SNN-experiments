library(ggplot2)
library(scales)

m_cmp_rslt <- read.table("m_cmp.csv", header = FALSE, sep = ",")
mtd_cmp_rslt <- read.table("mtd_cmp.csv", header = FALSE, sep = ",")
colnames(m_cmp_rslt) <- c("scenario", "m", "score", "method", "cov_kernel", "value")
colnames(mtd_cmp_rslt) <- c("scenario", "score", "method", "cov_kernel", "value")
mtd_cmp_rslt$m <- NA
mtd_cmp_rslt$m[mtd_cmp_rslt$method == "SNN"] <- 30
mtd_cmp_rslt$m[mtd_cmp_rslt$method == "SNN_order_maximin"] <- 30
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
all_rslt$submethod[ind_VT] <- "VMET[30]"
ind_TN <- which(all_rslt$method == "TN")
all_rslt$submethod[ind_TN] <- "MET"
ind_CB <- which(all_rslt$method == "CB")
all_rslt$submethod[ind_CB] <- "CSB"

m_vec <- sort(unique(m_cmp_rslt$m))
m_length <- length(m_vec)

# comparison between SNN and others ---------------------------------------
score <- "CRPS"
order <- NULL # "maximin" # scenario 1 and 2 do not have desc ordering results
if (is.null(order)) {
  mtd_vec <- c("CB", "SNN", "VT", "TN")
} else {
  mtd_vec <- c("CB", paste0("SNN", "_order_", order), "VT", "TN")
}
for (scenario in sort(unique(m_cmp_rslt$scenario))) {
  subset_mask <- all_rslt$score == score & all_rslt$scenario == scenario &
    all_rslt$method %in% mtd_vec
  subset_rslt <- all_rslt[subset_mask, , drop = FALSE]
  subset_rslt$submethod <- factor(subset_rslt$submethod, levels = c(paste0(
    "SNN",
    "[", m_vec, "]_known"
  ), "VMET[30]", "MET", "CSB", "SNN[30]_unknown"))
  ggplot(data = subset_rslt, mapping = aes(x = submethod, y = value)) +
    geom_boxplot(mapping = aes(fill = method), alpha = 0.5) +
    scale_y_continuous(name = score) +
    scale_x_discrete(labels = c(
      expression(SNN[10]),
      expression(SNN[20]),
      expression(SNN[30]),
      expression(SNN[40]),
      expression(SNN[50]),
      expression(VMET[30]),
      expression(MET),
      expression(CSB^"*"),
      expression(SNN[30]^"*")
    )) +
    theme(
      text = element_text(size = 16), legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
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
      "plots/performance_plot_scene_", scenario, "_", score,
      "_order_", order_ID, ".pdf"
    ),
    width = 5.5,
    height = 5
  )
}

# comparison of different orderings using SNN --------------------------------
score <- "CRPS"
mtd_vec <- c("SNN", "SNN_order_desc", "SNN_order_maximin")
scenario <- 3

subset_mask <- all_rslt$score == score & all_rslt$scenario == scenario &
  all_rslt$method %in% mtd_vec
subset_rslt <- all_rslt[subset_mask, , drop = FALSE]
subset_rslt$method <- factor(subset_rslt$method,
  levels = c("SNN", "SNN_order_desc", "SNN_order_maximin")
)
subset_rslt$m <- factor(subset_rslt$m, levels = sort(unique(subset_rslt$m)))
ggplot(data = subset_rslt, mapping = aes(x = m, y = value)) +
  geom_boxplot(mapping = aes(fill = method), alpha = 0.5) +
  scale_y_continuous(name = score) +
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
  paste0("plots/order_compare_scene_", scenario, "_", score, ".pdf"),
  width = 5,
  height = 5
)
