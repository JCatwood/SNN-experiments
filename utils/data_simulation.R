if (length(scene_ID) != 1L || is.na(scene_ID) || !scene_ID %in% 1:4) {
  stop("scene_ID must be one of 1, 2, 3, or 4.")
}

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
        row.names = FALSE, col.names = FALSE
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

# Scenario 2 is a mean-zero Gaussian process on a 100-by-100 grid with a
# nonstationary covariance constructed from a normalized random Gram matrix with
# a 0.1 nugget. 
# The GP field is censored below 1.
if (scene_ID == 2) {
  set.seed(123)
  tmp_vec <- seq(from = 0, to = 1, length.out = 100)
  locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
  mean_locs <- (locs[, 1] - 0.5)^2 + (locs[, 2] - 0.5)^2
  mat_tmp <- matrix(
    rnorm(100 * nrow(locs), mean = mean_locs, sd = sd(mean_locs)),
    nrow = 100, byrow = TRUE
  )
  covmat_tmp <- crossprod(mat_tmp)
  inv_sqrtdiag_covmat_tmp <- 1 / sqrt(diag(covmat_tmp))
  covmat <- outer(inv_sqrtdiag_covmat_tmp, inv_sqrtdiag_covmat_tmp) *
    covmat_tmp
  diag(covmat) <- diag(covmat) + 0.1
  rm(mean_locs, mat_tmp, covmat_tmp, inv_sqrtdiag_covmat_tmp)
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
                  row.names = FALSE, col.names = FALSE
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

# Scenario 3 is the mean-zero Matern 1.5 GP realized over [0, 1]^4.
# The GP field is censored below 1.
if (scene_ID == 3) {
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("Scenario 3 requires the lhs package.")
  }
  set.seed(123)
  locs <- lhs::randomLHS(n = 1e4, k = 4)
  cov_func <- GpGp::matern15_isotropic
  cov_parms <- c(1.0, 0.1, 0.0001)
  cov_name <- "matern15_isotropic"
  covmat <- cov_func(cov_parms, locs)
  N <- 20
  if (!file.exists("data/scenario_3")) {
    dir.create("data/scenario_3", recursive = TRUE)
  }
  if (all(file.exists(paste0("data/scenario_3/y", c(1:N), ".txt")))) {
    cat("Using previously generated GP realizations \n")
    y_list <- lapply(c(1:N), function(x) {
      as.vector(read.table(paste0("data/scenario_3/y", x, ".txt"), header = FALSE)[, 1])
    })
  } else {
    cat("Generating GP ...", "\n")
    L <- t(chol(covmat))
    y_list <- lapply(c(1:N), function(x) {
      as.vector(L %*% rnorm(nrow(locs)))
    })
    lapply(c(1:N), function(x) {
      write.table(y_list[[x]],
        file = paste0("data/scenario_3/y", x, ".txt"),
        row.names = FALSE, col.names = FALSE
      )
    })
    rm(L)
    cat("GP generated", "\n")
  }
  n <- nrow(locs)
  cens_lb <- rep(-Inf, n)
  cens_ub <- rep(1, n)
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
