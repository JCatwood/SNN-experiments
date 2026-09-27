library(R.utils)
offset <- 0.015

sample_func_VT_TN <- function(x1, x2, y1, y2, method = c("VT", "TN")) {
  mask_envelop <- locs[, 1] >= (x1 - offset) &
    locs[, 1] <= (x2 + offset) &
    locs[, 2] >= (y1 - offset) &
    locs[, 2] <= (y2 + offset)
  if (!any(mask_cens[mask_envelop])) return(list(ind = integer(), samp = matrix(numeric(), 0, n_samp)))
  locs_envelop <- locs[mask_envelop, , drop = FALSE]
  y_obs_envelop <- y_obs[mask_envelop]
  mask_cens_envelop <- mask_cens[mask_envelop]
  cens_ub_envelop <- cens_ub[mask_envelop]
  cens_lb_envelop <- cens_lb[mask_envelop]
  covmat_envelop <- covmat[mask_envelop, mask_envelop, drop = FALSE]
  locs_cens_envelop <- locs_envelop[mask_cens_envelop, , drop = FALSE]
  mask_inner <- locs_envelop[, 1] >= x1 &
    (locs_envelop[, 1] < x2 | x2 == 1) &
    locs_envelop[, 2] >= y1 & (locs_envelop[, 2] < y2 | y2 == 1)

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

  if (method[1] == "VT" && sum(mask_cens_envelop) > 1) {
    samp_envelop <- VeccTMVN::mvrandn(
      lower = cens_lb_envelop[mask_cens_envelop],
      upper = cens_ub_envelop[mask_cens_envelop],
      mean = cond_mean_cens_envelop,
      sigma = cond_covmat_cens_envelop,
      m = min(m, length(cond_mean_cens_envelop) - 1), N = n_samp
    )
  } else if (method[1] %in% c("TN", "VT")) {
    if (sum(mask_cens_envelop) > 2000) {
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
    (locs_cens_envelop[, 1] < x2 | x2 == 1) &
    locs_cens_envelop[, 2] >= y1 &
    (locs_cens_envelop[, 2] < y2 | y2 == 1)
  ind <- c(1:n)[mask_envelop][mask_inner & mask_cens_envelop]
  return(list(ind = ind, samp = samp_envelop[mask_cens_inner, , drop = FALSE]))
}

# Only dimension limits and timeouts cause subdivision; other errors surface.
sample_wrapper <- function(x1, x2, y1, y2, method = "VT", depth = 0) {
  inside <- locs[, 1] >= x1 - offset & locs[, 1] <= x2 + offset &
    locs[, 2] >= y1 - offset & locs[, 2] <= y2 + offset
  too_large <- sum(mask_cens[inside]) > 2000
  result <- if (too_large) NULL else tryCatch(
    R.utils::withTimeout(sample_func_VT_TN(x1, x2, y1, y2, method), timeout = 600),
    TimeoutException = function(e) NULL)
  if (!is.null(result)) return(result)
  if (depth == 4) stop("Partition limit reached; inspect this scenario before increasing it.")
  xs <- seq(x1, x2, length.out = 4)
  ys <- seq(y1, y2, length.out = 4)
  tiles <- expand.grid(x = 1:3, y = 1:3)
  result <- lapply(seq_len(nrow(tiles)), function(i) {
    a <- tiles$x[i]; b <- tiles$y[i]
    sample_wrapper(xs[a], xs[a + 1], ys[b], ys[b + 1], method, depth + 1)
  })
  list(ind = unlist(lapply(result, `[[`, "ind")),
    samp = do.call(rbind, lapply(result, `[[`, "samp")))
}
