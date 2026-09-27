source("../utils/local_moments.R")
args <- commandArgs(TRUE)
N <- if (length(args)) as.integer(args[1]) else 1000
reference <- length(args) > 1 && tolower(args[2]) == "true"
set.seed(123)
locs <- as.matrix(expand.grid(seq(0.025, 0.975, 0.05), seq(0.025, 0.975, 0.05)))
n <- nrow(locs)
Sigma <- fields::Matern(as.matrix(dist(locs)), range = 0.1, nu = 1.5)
y <- as.vector(t(chol(Sigma)) %*% rnorm(n))
cens <- y < 1
targets <- sample(which(cens), 10)
m_vec <- seq(10, 100, 10)
all <- local_moments(y, rep(-Inf, n), rep(1, n), Sigma, cens, targets,
  c(m_vec, if (reference) n), N, locs)
all$context <- "Observed and censored"
# u <- which(cens)
# only <- local_moments(y[u], rep(-Inf, length(u)), rep(1, length(u)),
#   Sigma[u, u], rep(TRUE, length(u)), match(targets, u),
#   c(m_vec, if (reference) length(u)), N, locs[u, ])
# only$index <- u[only$index]
# only$context <- "Censored only"
result <- rbind(all)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
write.csv(result, "results/local_moments_lowdim.csv", row.names = FALSE)
plot_moments(subset(result, m <= 100), "plots/local_moments_lowdim.pdf")
saveRDS(list(data = result, targets = targets, seed = 123, reference = reference,
  session = sessionInfo()), "results/local_moments_lowdim.rds")
