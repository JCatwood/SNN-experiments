# MET estimates of local moments; m counts the target and its neighbors.
local_moments <- function(y, lb, ub, Sigma, cens, targets, m_vec,
                          N = 1000, locs = NULL, seed = 123) {
  n <- length(y)
  draw <- function(j) {
    o <- j[!cens[j]]
    u <- j[cens[j]]
    mu <- rep(0, length(u))
    V <- Sigma[u, u, drop = FALSE]
    if (length(o)) {
      B <- solve(Sigma[o, o, drop = FALSE], Sigma[o, u, drop = FALSE])
      mu <- as.vector(crossprod(B, y[o]))
      V <- V - crossprod(B, Sigma[o, u, drop = FALSE])
    }
    X <- TruncatedNormal::rtmvnorm(N, mu, (V + t(V)) / 2, lb[u], ub[u])
    colnames(X) <- u
    X
  }
  out <- list()
  for (m in unique(pmin(m_vec, n))) {
    set.seed(seed + m)
    full <- if (m == n) draw(seq_len(n)) else NULL
    for (i in targets) {
      distance <- if (is.null(locs)) {
        1 - abs(Sigma[i, ] / sqrt(Sigma[i, i] * diag(Sigma)))
      } else rowSums(sweep(locs, 2, locs[i, ], "-")^2)
      j <- c(i, setdiff(order(distance), i))[seq_len(m)]
      x <- if (m == n) full[, as.character(i)] else draw(j)[, as.character(i)]
      v <- var(x)
      se_v <- sqrt(max(0, mean((x - mean(x))^4) - v^2) / N)
      out[[length(out) + 1]] <- data.frame(index = i, m = m, N = N,
        mean = mean(x), variance = v, mean_se = sqrt(v / N), variance_se = se_v)
    }
  }
  do.call(rbind, out)
}

plot_moments <- function(x, file) {
  long <- rbind(transform(x, moment = "Mean", value = mean, se = mean_se),
                transform(x, moment = "Variance", value = variance, se = variance_se))
  p <- ggplot2::ggplot(long, ggplot2::aes(m, value, color = factor(index))) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = value - 1.96 * se,
      ymax = value + 1.96 * se, fill = factor(index)), alpha = 0.2, color = NA) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 0.8) +
    ggplot2::facet_wrap(~ moment, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "Neighborhood size m", y = NULL, color = "Index", fill = "Index") +
    ggplot2::theme_bw(base_size = 11) + ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file, p, width = 7, height = 3)
}
