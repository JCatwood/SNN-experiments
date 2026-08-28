rm(list = ls())

# command-line arguments ---------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript highdim_truth_unknown.R [k scene_ID run_CB]",
  "  k: positive integer indexing the GP realization",
  "  scene_ID: positive integer identifying the simulation scenario",
  "  run_CB: true or false (default: true)",
  sep = "\n"
)

if (any(args %in% c("-h", "--help"))) {
  cat(usage, "\n")
  quit(save = "no", status = 0)
}

parse_bool <- function(x) {
  value <- tolower(x)
  if (value %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop("run_CB must be true or false.\n", usage, call. = FALSE)
}

k <- 1
scene_ID <- 1
run_CB <- TRUE
if (length(args) > 0) {
  if (!length(args) %in% c(2, 3)) {
    stop(usage, call. = FALSE)
  }
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2]) # simulation scenario ID
  if (length(args) == 3) {
    run_CB <- parse_bool(args[3])
  }
}

if (anyNA(c(k, scene_ID)) || any(c(k, scene_ID) < 1)) {
  stop("k and scene_ID must be positive integers.\n", usage, call. = FALSE)
}
if (run_CB && !scene_ID %in% c(1, 2)) {
  stop(
    "The CSB comparison is restricted to Scenarios 1 and 2.",
    call. = FALSE
  )
}

library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(scoringRules)

if (run_CB && !requireNamespace("CensSpBayes", quietly = TRUE)) {
  stop(
    "run_CB=true requires the optional CensSpBayes package. ",
    "Run with run_CB=false to skip this experiment.",
    call. = FALSE
  )
}

# simulation settings ------------------------------------------------------
set.seed(123)
m <- 30 # number of nearest neighbors
n_samp <- 50 # samples generated for posterior inference
# CensSpBayes
n_burn <- 20000
n_iter_MC <- 25000
thin <- 5

# data simulation ----------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
y_test <- y_test_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- cens_ub[mask_cens] # CensSpBayes does not allow NA
if (!exists("cov_name")) {
  cov_name <- "matern15_isotropic"
}
L <- t(chol(covmat))

source("../utils/score_output.R")

# CensSpBayes ------------------------------
if (run_CB) {
  bgn_time <- Sys.time()
  inla.mats <- CensSpBayes::create_inla_mats(
    S = locs, S.pred = locs[mask_cens, ],
    offset = c(0.01, 0.2),
    cutoff = 0.05,
    max.edge = c(0.01, 0.1)
  )
  X.obs <- matrix(1, nrow(locs), 1)
  X.pred <- matrix(1, sum(mask_cens), 1)
  cat("CB sampling begins...\n")
  ret_obj <- CensSpBayes::CensSpBayes(
    Y = y_obs, S = locs, X = X.obs,
    cutoff.Y = cens_ub,
    S.pred = locs[mask_cens, ], X.pred = X.pred,
    inla.mats = inla.mats,
    rho.init = 0.1, rho.upper = 5,
    iters = n_iter_MC, burn = n_burn, thin = thin, ret_samp = TRUE
  )
  y_samp_CB <- matrix(y_obs,
    nrow = length(y_obs),
    ncol = ncol(ret_obj$Y.pred.samp), byrow = FALSE
  )
  y_samp_CB[mask_cens, ] <- ret_obj$Y.pred.samp
  cat("CB sampling done\n")
  end_time <- Sys.time()
  time_CB <- difftime(end_time, bgn_time, units = "secs")[[1]]

  kriging_score_output(y_samp_CB, y_test, time_CB,
    scene_ID = scene_ID, method = "CB", parms = "unknown"
  )
}

# save data for heatmap -------------------------
if (k == 1) {
  if (!file.exists("samples")) {
    dir.create("samples")
  }
  if (run_CB) {
    y_all <- rep(NA, length(ind_train) + length(ind_test))
    y_all[ind_train] <- y_samp_CB[, 1]
    y_all[ind_test] <- y_test
    write.table(y_all,
      file = paste0(
        "samples/samp_cmp_known_CB_scene",
        scene_ID, ".csv"
      ),
      row.names = F, col.names = F
    )
  }
}
