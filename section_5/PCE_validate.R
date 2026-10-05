# Run from section_5: Rscript PCE_validate.R TX 123 1000 SNN,VMET,MET,LOD
source("../utils/predict_marginals.R")
args <- commandArgs(TRUE)
region <- if (length(args)) toupper(args[1]) else "TX"
seed <- if (length(args) > 1) as.integer(args[2]) else 123
N <- if (length(args) > 2) as.integer(args[3]) else 50
methods <- if (length(args) > 3) strsplit(args[4], ",")[[1]] else
  if (region == "TX") c("SNN", "VMET", "MET", "LOD") else c("SNN", "LOD")
stopifnot(region %in% c("TX", "US"), all(methods %in% c("SNN", "VMET", "MET", "LOD", "CSB")))
load("PCE.RData")
dat <- data.PCE.censored
region_mask <- rep(TRUE, nrow(dat))
if (region == "TX") {
  states <- spData::us_states
  pts <- sf::st_as_sf(dat, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  texas <- states[states$NAME == "Texas", ]
  in_tx <- lengths(sf::st_intersects(pts, sf::st_transform(texas, 4326))) > 0
  # keep <- in_tx | !dat$left_censored
  keep <- in_tx
  region_mask <- in_tx[keep]
  dat <- dat[keep, ]
}
y <- log(dat$result_va + 1e-8)
b <- log(dat$detection_level + 1e-8)
cens <- as.logical(dat$left_censored)
x <- cbind(dat$lon, dat$lat, as.numeric(dat$startDate))
# Hold out complete locations, including repeated visits to the same well.
group <- paste(dat$lon, dat$lat, sep = ":")
set.seed(seed)
held <- sample(unique(group), max(1, round(0.2 * length(unique(group)))))
test <- which(group %in% held & region_mask)
train <- which(!group %in% held)
center <- mean(y[train][!cens[train]])
scale_y <- sd(y[train][!cens[train]])
y <- (y - center) / scale_y
b <- (b - center) / scale_y
origin <- apply(x[train, ], 2, min)
span <- apply(x[train, ], 2, function(z) diff(range(z)))
x <- sweep(sweep(x, 2, origin, "-"), 2, pmax(span, 1), "/")
yt <- y[train]
yt[cens[train]] <- NA
bt <- b[train]
ct <- cens[train]
xt <- x[train, , drop = FALSE]
if (region == "TX") {
  fit <- seq_len(length(train))  
} else {
  fit <- GpGp::order_maxmin(xt)[seq_len(max(51, round(0.2 * length(train))))]
  fit <- fit[!is.na(fit)]
}
nll <- function(lp) {
  set.seed(123)
  -VeccTMVN::loglk_censor_MVN(xt[fit, ], which(ct[fit]), yt[fit], bt[fit],
    "matern15_scaledim", exp(lp), m = min(50, length(fit) - 1), NLevel2 = 1000)
}
fit_time <- system.time(opt <- optim(log(c(1, .01, .01, .01, .005)), nll,
  control = list(maxit = 200)))[[3]]
# if (opt$convergence != 0) stop("Training-only covariance fit did not converge; inspect optimizer.")
par <- exp(opt$par)
# Range-scaled geometry for both SNN and the common marginal prediction step.
xs <- sweep(x, 2, par[2:4], "/")
ps <- par
ps[2:4] <- 1
u <- which(ct)
o <- which(!ct)
if (any(methods %in% c("VMET", "MET"))) {
  K <- GpGp::matern15_scaledim(par, xt)
  B <- solve(K[o, o], K[o, u, drop = FALSE])
  conditional_mean <- as.vector(crossprod(B, yt[o]))
  V <- K[u, u, drop = FALSE] - crossprod(B, K[o, u, drop = FALSE])
  V <- (V + t(V)) / 2
}
dir.create("results", showWarnings = FALSE)
records <- list()
for (method in methods) {
  started <- proc.time()[[3]]
  set.seed(seed)
  draws <- matrix(yt, length(train), N)
  if (method == "SNN") {
    draws <- vapply(seq_len(N), function(s) nntmvn::rptmvn(yt,
      rep(-Inf, length(train)), bt, ct, m = 50, locs = xs[train, ],
      cov_name = "matern15_scaledim", cov_parm = ps, ordering = 2,
      seed = seed + s), numeric(length(train)))
  } else if (method == "MET") {
    draws[u, ] <- t(TruncatedNormal::rtmvnorm(N, conditional_mean, V,
      rep(-Inf, length(u)), bt[u]))
  } else if (method == "VMET") {
    draws[u, ] <- VeccTMVN::mvrandn(lower = rep(-Inf, length(u)), upper = bt[u],
      mean = conditional_mean, sigma = V, N = N, m = min(50, length(u) - 1))
  } else if (method == "LOD") {
    draws[u, ] <- bt[u]
  }
  if (method == "CSB") {
    # CSB fits its own spatial model; test outcomes/limits are never conditioning data.
    S <- x[train, 1:2, drop = FALSE]
    Sp <- x[test, 1:2, drop = FALSE]
    mats <- CensSpBayes::create_inla_mats(S = S, S.pred = Sp,
      offset = c(.01, .2), cutoff = .05, max.edge = c(.01, .1))
    obs <- yt
    obs[ct] <- bt[ct]
    cb <- CensSpBayes::CensSpBayes(Y = obs, S = S, X = matrix(1, nrow(S), 1),
      cutoff.Y = bt, S.pred = Sp, X.pred = matrix(1, nrow(Sp), 1),
      inla.mats = mats, rho.init = .1, rho.upper = 5,
      iters = 20000 + 5 * 1000, burn = 20000, thin = 5, ret_samp = TRUE)
    pred <- cb$Y.pred.samp
    mean_pred <- rowMeans(pred)
    probability <- rowMeans(pred <= b[test])
  } else {
    p <- predict_marginals(xs[train, ], xs[test, ], draws, ps, seed = seed)
    pred <- p$draws
    mean_pred <- rowMeans(p$mean)
    probability <- rowMeans(pnorm((b[test] - p$mean) / sqrt(p$variance)))
  }
  elapsed <- proc.time()[[3]] - started
  observed <- !cens[test] & is.finite(y[test])
  eligible <- is.finite(b[test])
  raw <- pred * scale_y + center
  actual <- y[test] * scale_y + center
  interval <- t(apply(raw, 1, quantile, c(.05, .95)))
  # metrics <- c(
  #   RMSE = sqrt(mean((actual[observed] -
  #     (mean_pred[observed] * scale_y + center))^2)),
  #   CRPS = mean(scoringRules::crps_sample(actual[observed], raw[observed, , drop = FALSE])),
  #   Brier = mean((probability[eligible] - cens[test][eligible])^2),
  #   Coverage90 = mean(actual[observed] >= interval[observed, 1] &
  #     actual[observed] <= interval[observed, 2]))
  metrics <- c(Brier = mean((probability[eligible] - cens[test][eligible])^2))
  records[[method]] <- data.frame(region, seed, method, metric = names(metrics),
    value = as.numeric(metrics), n_observed = sum(observed), n_status = sum(eligible),
    n_status_observed = sum(eligible & !cens[test]), n_status_censored = sum(eligible & cens[test]),
    fit_seconds = if (method == "CSB") 0 else fit_time, prediction_seconds = elapsed)
  saveRDS(list(region = region, seed = seed, method = method, train = train,
    test = test, rows = rownames(dat), groups = group, center = center, scale = scale_y, optimizer = opt,
    predictions = raw, probability_censored = probability, observed = observed,
    censored = cens[test], upper = b[test] * scale_y + center,
    session = sessionInfo()), paste0("results/PCE_validation_", region, "_", seed, "_", method, ".rds"))
}
write.csv(do.call(rbind, records), paste0("results/PCE_validation_", region, "_", seed, ".csv"),
  row.names = FALSE)
