library(R.utils)
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
        timeout = 600
      )
    },
    TimeoutException = function(foo) {
      cat(
        method, "sampling in [", x1, x2, "] X [", y1, y2, "] did not finish",
        "within 10 minutes\n"
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
