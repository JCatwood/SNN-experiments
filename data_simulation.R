# Scenario 1 is the mean-zero Matern 1.5 GP realized over [0, 1]^2.
# The GP field is censored below 1.
if (scene_ID == 1) {
  set.seed(123)
  tmp_vec <- seq(from = 0, to = 1, length.out = 100)
  locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
  cov_func <- GpGp::matern15_isotropic
  cov_parms <- c(1.0, 0.03, 0.0001)
  cov_name <- "matern15_isotropic"
  covmat <- cov_func(cov_parms, locs)
  N <- 20
  if (!file.exists("data/scenario_1")) {
    dir.create("data/scenario_1", recursive = TRUE)
  }
  if (all(file.exists(paste0("data/scenario_1/y", c(1:N), ".txt")))) {
    cat("Using previously generated GP realizations \n")
    y_list <- lapply(c(1:N), function(x) {
      as.vector(read.table(paste0("data/scenario_1/y", x, ".txt"), header = FALSE)[, 1])
    })
  } else {
    cat("Generating GP ...", "\n")
    L <- t(chol(covmat))
    y_list <- lapply(c(1:N), function(x) {
      as.vector(L %*% rnorm(nrow(locs)))
    })
    lapply(c(1:N), function(x) {
      write.table(y_list[[x]],
        file = paste0("data/scenario_1/y", x, ".txt"),
        row.names = FALSE, col.names = FALSE,
      )
    })
    rm(L)
    cat("GP generated", "\n")
  }
  n <- nrow(locs)
  cens_lb <- rep(-Inf, n)
  cens_ub <- rep(1, n)
  rm(tmp_vec)
}

# Scenario 2 is the Himmelblau's function (https://en.wikipedia.org/wiki/Himmelblau%27s_function)
# defined in the region of [-2.5, 2.5]^2.
# The values from the Himmelblau's function is normalized to have standard
# deviation of 1, after which normal noise with standard deviation 0.1 is added.
# The random field is censored above 0.
if (scene_ID == 2) {
  set.seed(123)
  tmp_vec <- seq(from = 0, to = 1, length.out = 100)
  locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
  y <- (((locs[, 1] - 0.5) * 5)^2 + ((locs[, 2] - 0.5) * 5) - 11) +
    (((locs[, 1] - 0.5) * 10) + ((locs[, 2] - 0.5) * 10)^2 - 7)
  y <- y / sd(y) + rnorm(length(y), sd = 0.1)
  n <- nrow(locs)
  cens_lb <- rep(0, n)
  cens_ub <- rep(Inf, n)
  rm(tmp_vec)
}
