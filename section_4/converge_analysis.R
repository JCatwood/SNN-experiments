source("../utils/local_moments.R")
args <- commandArgs(TRUE)
k <- if (length(args)) as.integer(args[1]) else 1
scene_ID <- if (length(args) > 1) as.integer(args[2]) else 1
N <- if (length(args) > 2) as.integer(args[3]) else 1000
source("../utils/data_simulation.R")
y <- y_list[[k]]
cens <- y > cens_lb & y < cens_ub
set.seed(123)
targets <- sample(which(cens), 10)
# Correlation neighborhoods for periodic and random covariances.
metric_locs <- if (scene_ID %in% c(2)) NULL else locs
result <- local_moments(y, cens_lb, cens_ub, covmat, cens, targets,
  seq(10, 100, 10), N, metric_locs)
result$context <- paste("Scenario", scene_ID)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
stem <- paste0("local_moments_scene_", scene_ID, "_seed_", k)
write.csv(result, paste0("results/", stem, ".csv"), row.names = FALSE)
plot_moments(result, paste0("plots/", stem, ".pdf"))
saveRDS(list(data = result, targets = targets, session = sessionInfo()),
  paste0("results/", stem, ".rds"))
