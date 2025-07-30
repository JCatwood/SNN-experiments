library(doParallel)
library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(scoringRules)

# simulation settings ------------------------------
rm(list = ls())
set.seed(123)
m <- 30 # number of nearest neighbors
n_samp <- 50 # samples generated for posterior inference
run_SNN <- TRUE
run_VT <- TRUE
run_TN <- TRUE
use_parallel <- FALSE
plot_heatmap <- FALSE
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2])
  use_snn_order <- as.integer(args[3]) # 0, 1, 2
} else {
  k <- 1
  scene_ID <- 1
  use_snn_order <- 0
}

# data simulation ----------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)
y_obs <- y
y_obs[mask_cens] <- NA

# nntmvn ---------------------------------------
if (run_SNN) {
  order <- 1:n
  y_obs_order <- y_obs[order]
  cens_ub_order <- cens_ub[order]
  cens_lb_order <- cens_lb[order]
  mask_cens_order <- mask_cens[order]
  locs_order <- locs[order, , drop = FALSE]
  covmat_order <- covmat[order, order]
  cat("SNN sampling...\n")
  bgn_time <- Sys.time()
  if (use_parallel) {
    ncores <- 4
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    y_samp_SNN_order <- foreach(i = 1:n_samp, .packages = c("nntmvn")) %dopar% {
      nntmvn::rptmvn(y_obs_order, cens_lb_order, cens_ub_order,
        mask_cens_order,
        m = m,
        covmat = covmat_order, locs = locs_order, ordering = use_snn_order,
        seed = i
      )
    }
    stopCluster(cl)
  } else {
    y_samp_SNN_order <- lapply(1:n_samp, function(seed_id) {
      nntmvn::rptmvn(y_obs_order, cens_lb_order, cens_ub_order,
        mask_cens_order,
        m = m,
        covmat = covmat_order, locs = locs_order, ordering = use_snn_order,
        seed = seed_id
      )
    })
  }
  end_time <- Sys.time()
  cat("SNN sampling is done\n")
  time_SNN <- difftime(end_time, bgn_time, units = "secs")[[1]]
  rev_order <- 1:n
  rev_order[order] <- 1:n
  y_samp_SNN <- matrix(unlist(y_samp_SNN_order), n,
    n_samp,
    byrow = FALSE
  )[rev_order, , drop = FALSE]
  y_pred_SNN <- rowMeans(y_samp_SNN)
  y_pred_cens_SNN <- y_pred_SNN[mask_cens]
  cat(
    "> ", scene_ID, ", RMSE, SNN, known, ",
    sqrt(mean((y[mask_cens] - y_pred_cens_SNN)^2)), "\n"
  )
  cat(
    "> ", scene_ID, ", CRPS, SNN, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens],
      dat = y_samp_SNN[mask_cens, , drop = FALSE]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ", time, SNN, known, ",
    time_SNN, "\n"
  )
}

# VeccTMVN and TruncatedNormal -----------------------
offset <- 0.015
sample_func_VT_TN <- function(x1, x2, y1, y2, method = c("VT", "TN")) {
  mask_envelop <- locs[, 1] >= (x1 - offset) &
    locs[, 1] <= (x2 + offset) &
    locs[, 2] >= (y1 - offset) &
    locs[, 2] <= (y2 + offset)
  locs_envelop <- locs[mask_envelop, , drop = FALSE]
  y_obs_envelop <- y_obs[mask_envelop]
  mask_cens_envelop <- mask_cens[mask_envelop]
  cens_ub_envelop <- cens_ub[mask_envelop]
  cens_lb_envelop <- cens_lb[mask_envelop]
  covmat_envelop <- covmat[mask_envelop, mask_envelop, drop = FALSE]
  locs_cens_envelop <- locs_envelop[mask_cens_envelop, , drop = FALSE]
  mask_inner <- locs_envelop[, 1] >= x1 &
    locs_envelop[, 1] <= x2 &
    locs_envelop[, 2] >= y1 & locs_envelop[, 2] <= y2

  cat(
    "Dimension of the TMVN distribution to be sampled from is",
    sum(mask_cens_envelop), "\n"
  )

  if (sum(mask_cens_envelop) == length(y_obs_envelop)) {
    cond_mean_cens_envelop <- rep(0, sum(mask_envelop))
    cond_covmat_cens_envelop <- covmat_envelop
  } else {
    tmp_mat <- solve(
      covmat_envelop[!mask_cens_envelop, !mask_cens_envelop, drop = FALSE],
      covmat_envelop[!mask_cens_envelop, mask_cens_envelop, drop = FALSE]
    )
    cond_covmat_cens_envelop <-
      covmat_envelop[mask_cens_envelop, mask_cens_envelop, drop = FALSE] -
      covmat_envelop[mask_cens_envelop, !mask_cens_envelop, drop = FALSE] %*%
      tmp_mat
    cond_mean_cens_envelop <- as.vector(
      t(y_obs_envelop[!mask_cens_envelop]) %*% tmp_mat
    )
  }
  # guarantee symmetry
  cond_covmat_cens_envelop[lower.tri(cond_covmat_cens_envelop)] <-
    t(cond_covmat_cens_envelop)[lower.tri(cond_covmat_cens_envelop)]

  if (method[1] == "VT") {
    samp_envelop <- VeccTMVN::mvrandn(
      lower = cens_lb_envelop[mask_cens_envelop],
      upper = cens_ub_envelop[mask_cens_envelop],
      mean = cond_mean_cens_envelop,
      sigma = cond_covmat_cens_envelop,
      m = min(m, length(cond_mean_cens_envelop) - 1), N = n_samp
    )
  } else if (method[1] == "TN") {
    if (sum(mask_cens_envelop) > 1500) {
      stop("Input dimension for TN is too high\n")
    }
    samp_envelop <- t(TruncatedNormal::rtmvnorm(
      n_samp, cond_mean_cens_envelop,
      cond_covmat_cens_envelop, cens_lb_envelop[mask_cens_envelop],
      cens_ub_envelop[mask_cens_envelop]
    ))
  } else {
    stop("invalid method option\n")
  }


  locs_cens_envelop <- locs_envelop[mask_cens_envelop, , drop = FALSE]
  mask_cens_inner <- locs_cens_envelop[, 1] >= x1 &
    locs_cens_envelop[, 1] <= x2 &
    locs_cens_envelop[, 2] >= y1 &
    locs_cens_envelop[, 2] <= y2
  ind <- c(1:n)[mask_envelop][mask_inner & mask_cens_envelop]
  return(list(ind = ind, samp = samp_envelop[mask_cens_inner, ]))
}

library(R.utils)
sample_wrapper <- function(x1, x2, y1, y2, method = c("VT", "TN")) {
  cat(method, "sampling [", x1, x2, "] X [", y1, y2, "]...\n")
  ret_obj <- tryCatch(
    {
      R.utils::withTimeout(
        {
          ret_obj <- sample_func_VT_TN(x1, x2, y1, y2, method)
          cat("Done", method, "sampling [", x1, x2, "] X [", y1, y2, "]\n")
          ret_obj
        },
        timeout = 300
      )
    },
    TimeoutException = function(foo) {
      cat(
        method, "sampling in [", x1, x2, "] X [", y1, y2, "] did not finish",
        "within 5 minutes\n"
      )
      x_span <- (x2 - x1) / 3
      y_span <- (y2 - y1) / 3
      ret_obj1 <- sample_wrapper(x1, x1 + x_span, y1, y1 + y_span, method)
      ret_obj2 <- sample_wrapper(x1, x1 + x_span, y1 + y_span, y1 + 2 * y_span, method)
      ret_obj3 <- sample_wrapper(x1, x1 + x_span, y1 + 2 * y_span, y2, method)
      ret_obj4 <- sample_wrapper(x1 + x_span, x1 + 2 * x_span, y1, y1 + y_span, method)
      ret_obj5 <- sample_wrapper(x1 + x_span, x1 + 2 * x_span, y1 + y_span, y1 + 2 * y_span, method)
      ret_obj6 <- sample_wrapper(x1 + x_span, x1 + 2 * x_span, y1 + 2 * y_span, y2, method)
      ret_obj7 <- sample_wrapper(x1 + 2 * x_span, x2, y1, y1 + y_span, method)
      ret_obj8 <- sample_wrapper(x1 + 2 * x_span, x2, y1 + y_span, y1 + 2 * y_span, method)
      ret_obj9 <- sample_wrapper(x1 + 2 * x_span, x2, y1 + 2 * y_span, y2, method)

      ret_obj <- list(
        ind = c(
          ret_obj1$ind, ret_obj2$ind, ret_obj3$ind,
          ret_obj4$ind, ret_obj5$ind, ret_obj6$ind,
          ret_obj7$ind, ret_obj8$ind, ret_obj9$ind
        ),
        samp = rbind(
          ret_obj1$samp, ret_obj2$samp, ret_obj3$samp,
          ret_obj4$samp, ret_obj5$samp, ret_obj6$samp,
          ret_obj7$samp, ret_obj8$samp, ret_obj9$samp
        )
      )
      ret_obj
    },
    error = function(e) {
      message("Caught error: ", e$message)
      x_span <- (x2 - x1) / 3
      y_span <- (y2 - y1) / 3
      ret_obj1 <- sample_wrapper(x1, x1 + x_span, y1, y1 + y_span, method)
      ret_obj2 <- sample_wrapper(x1, x1 + x_span, y1 + y_span, y1 + 2 * y_span, method)
      ret_obj3 <- sample_wrapper(x1, x1 + x_span, y1 + 2 * y_span, y2, method)
      ret_obj4 <- sample_wrapper(x1 + x_span, x1 + 2 * x_span, y1, y1 + y_span, method)
      ret_obj5 <- sample_wrapper(x1 + x_span, x1 + 2 * x_span, y1 + y_span, y1 + 2 * y_span, method)
      ret_obj6 <- sample_wrapper(x1 + x_span, x1 + 2 * x_span, y1 + 2 * y_span, y2, method)
      ret_obj7 <- sample_wrapper(x1 + 2 * x_span, x2, y1, y1 + y_span, method)
      ret_obj8 <- sample_wrapper(x1 + 2 * x_span, x2, y1 + y_span, y1 + 2 * y_span, method)
      ret_obj9 <- sample_wrapper(x1 + 2 * x_span, x2, y1 + 2 * y_span, y2, method)

      ret_obj <- list(
        ind = c(
          ret_obj1$ind, ret_obj2$ind, ret_obj3$ind,
          ret_obj4$ind, ret_obj5$ind, ret_obj6$ind,
          ret_obj7$ind, ret_obj8$ind, ret_obj9$ind
        ),
        samp = rbind(
          ret_obj1$samp, ret_obj2$samp, ret_obj3$samp,
          ret_obj4$samp, ret_obj5$samp, ret_obj6$samp,
          ret_obj7$samp, ret_obj8$samp, ret_obj9$samp
        )
      )
      ret_obj
    }
  )
  return(ret_obj)
}

if (run_VT) {
  y_samp_VT <- matrix(y_obs,
    nrow = length(y_obs),
    ncol = n_samp, byrow = FALSE
  )
  bgn_time <- Sys.time()
  set.seed(123)
  ret_obj <- sample_wrapper(0, 1, 0, 1, "VT")
  y_samp_VT[ret_obj$ind, ] <- ret_obj$samp
  end_time <- Sys.time()
  time_VT <- difftime(end_time, bgn_time, units = "secs")[[1]]
  y_pred_VT <- rowMeans(y_samp_VT)
  cat(
    "> ", scene_ID, ", RMSE, VT, known, ",
    sqrt(mean((y[mask_cens] - y_pred_VT[mask_cens])^2)), "\n"
  )
  cat(
    "> ", scene_ID, ", CRPS, VT, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens], dat = y_samp_VT[mask_cens, , drop = FALSE]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ", time, VT, known, ",
    time_VT, "\n"
  )
}

if (run_TN) {
  y_samp_TN <- matrix(y_obs,
    nrow = length(y_obs),
    ncol = n_samp, byrow = FALSE
  )
  bgn_time <- Sys.time()
  set.seed(123)
  ret_obj <- sample_wrapper(0, 1, 0, 1, "TN")
  y_samp_TN[ret_obj$ind, ] <- ret_obj$samp
  end_time <- Sys.time()
  time_TN <- difftime(end_time, bgn_time, units = "secs")[[1]]
  y_pred_TN <- rowMeans(y_samp_TN)
  cat(
    "> ", scene_ID, ", RMSE, TN, known, ",
    sqrt(mean((y[mask_cens] - y_pred_TN[mask_cens])^2)), "\n"
  )
  cat(
    "> ", scene_ID, ", CRPS, TN, known, ",
    mean(scoringRules::crps_sample(
      y = y[mask_cens], dat = y_samp_TN[mask_cens, , drop = FALSE]
    )), "\n"
  )
  cat(
    "> ", scene_ID, ", time, TN, known, ",
    time_TN, "\n"
  )
}

# heatmap -------------------------
if (plot_heatmap) {
  library(fields)
  zlim <- range(y_samp_SNN, y_samp_TN, y_samp_VT)
  if (!file.exists("plots")) {
    dir.create("plots")
  }
  pdf(
    file = paste0("plots/sim_data_scene", scene_ID, "_true.pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y, sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0(
      "plots/samp_cmp_known_SNN_scene", scene_ID, "_m", m,
      "_order", use_snn_order, ".pdf"
    ),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_SNN[, 1], sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0("plots/samp_cmp_known_VT_scene", scene_ID, ".pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_VT[, 1], sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
  pdf(
    file = paste0("plots/samp_cmp_known_TN_scene", scene_ID, ".pdf"),
    width = 5, height = 5
  )
  fields::image.plot(matrix(y_samp_TN[, 1], sqrt(n), sqrt(n)), zlim = zlim)
  dev.off()
}
