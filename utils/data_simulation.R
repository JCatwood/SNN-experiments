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

# Scenario 2 is a 2d periodic covariance function over [0, 1]^2.
# The GP field is censored below 1.
if (scene_ID == 2) {
  set.seed(123)
  tmp_vec <- seq(from = 0, to = 1, length.out = 100)
  locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
  # define a 2d periodic covariance function
  cov_func <- function(cov_parms, locs) {
    variance <- cov_parms[1]
    period <- cov_parms[2]
    range <- cov_parms[3]
    nugget <- cov_parms[4]
    dist_mat_x <- as.matrix(dist(locs[, 1, drop = F]))
    dist_mat_y <- as.matrix(dist(locs[, 2, drop = F]))
    variance * exp(-2 * sin(pi * dist_mat_x / period)^2 / range^2) *
      exp(-2 * sin(pi * dist_mat_y / period)^2 / range^2) +
      diag(rep(nugget * variance, nrow(locs)))
  }
  cov_parms <- c(1.0, 0.5, 0.3, 0.0001)
  # cov_name used for model fitting, not the true covariance model
  cov_name <- "matern15_isotropic"
  covmat <- cov_func(cov_parms, locs)
  N <- 20
  if (!file.exists("data/scenario_2")) {
    dir.create("data/scenario_2", recursive = TRUE)
  }
  if (all(file.exists(paste0("data/scenario_2/y", c(1:N), ".txt")))) {
    cat("Using previously generated GP realizations \n")
    y_list <- lapply(c(1:N), function(x) {
      as.vector(read.table(paste0("data/scenario_2/y", x, ".txt"), header = FALSE)[, 1])
    })
  } else {
    cat("Generating GP ...", "\n")
    L <- t(chol(covmat))
    y_list <- lapply(c(1:N), function(x) {
      as.vector(L %*% rnorm(nrow(locs)))
    })
    lapply(c(1:N), function(x) {
      write.table(y_list[[x]],
        file = paste0("data/scenario_2/y", x, ".txt"),
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

# Scenario 3 is the same as Scenario 1 but with ub generated from N(0, 1)
# The GP field is censored below 1.
if (scene_ID == 3) {
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
  cens_ub <- rnorm(n)
  rm(tmp_vec)
}

# split data into training and testing -----------------------
set.seed(123)
n_test <- round(n * 0.2)
ind_test <- sample(1 : n, n_test, replace = FALSE)
ind_train <- c(1 : n)[-ind_test]
covmat_train_test <- covmat[ind_train, ind_test, drop = FALSE]
covmat_test <- covmat[ind_test, ind_test,  drop = FALSE]
covmat <- covmat[ind_train, ind_train, drop = FALSE]
locs_test <- locs[ind_test, , drop = FALSE]
locs <- locs[ind_train, , drop = FALSE]
y_test_list <- lapply(y_list, function(x) { x[ind_test] })
y_list <- lapply(y_list, function(x) { x[ind_train] })
cens_lb <- cens_lb[ind_train]
cens_ub <- cens_ub[ind_train]
n <- n - n_test