library(dplyr)
library(tidyr)
files <- list.files("results/", "\\.csv$", full.names = TRUE)
if (!length(files)) stop("Run the revised experiments first.")
if (!dir.exists("results/output_tables")) {
  dir.create("results/output_tables", recursive = TRUE)
}
x <- bind_rows(lapply(files, read.csv))
stopifnot(all(x$scoring_version == 2))
s <- x %>% group_by(scenario, method, cov_kernel, m, score) %>%
  summarise(mean = mean(value), se = sd(value) / sqrt(n()), replicates = n(), .groups = "drop")
table <- s %>% 
  mutate(mean = if_else(method %in% c("MET", "VT") & score == "time", mean + 600, mean)) %>%
  # mutate(cell = sprintf("%.4g (SE %.2g; n=%d)", mean, se, replicates)) %>%
  mutate(cell = sprintf("%.4g (%.2g)", mean, se)) %>%
  select(scenario, method, cov_kernel, m, score, cell) %>%
  pivot_wider(names_from = score, values_from = cell) %>%
  filter(method %in% c("CSB", "MET", "VT", "SNN_order_maximin")) %>%
  mutate(method = ifelse(method == "SNN_order_maximin", "SNN (maximin)", method)) %>%
  mutate(method = ifelse(method == "VT", "VMET", method)) %>%
  rename(Covariance = cov_kernel) %>%
  mutate(Covariance = ifelse(Covariance == "known", "true", "estimated"))

# md table  
header <- paste0("| ", paste(names(table), collapse = " | "), " |")
rows <- apply(table, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
writeLines(c(header, paste0("|", paste(rep("---", ncol(table)), collapse = "|"), "|"), rows),
  "results/output_tables/performance_table.md")

# tex table
header <- paste0(
  paste(names(table), collapse = " & "),
  " \\\\"
)
rows <- apply(table, 1, function(z) {
  paste0(paste(z, collapse = " & "), " \\\\")
})
latex_table <- c(
  "\\begin{tabular}{",
  paste(rep("c", ncol(table)), collapse = ""),
  "}",
  "\\hline",
  header,
  "\\hline",
  rows,
  "\\hline",
  "\\end{tabular}"
)
writeLines(
  latex_table,
  "results/output_tables/performance_table.tex"
)