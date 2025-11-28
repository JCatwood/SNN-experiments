library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(ggplot2)
library(nntmvn)

# lowdim example ------------------------------
rm(list = ls())
set.seed(123)
tmp_vec <- seq(from = 0.025, to = 0.975, by = 0.05)
locs <- as.matrix(expand.grid(tmp_vec, tmp_vec))
rm(tmp_vec)
n <- nrow(locs)
m <- 30
range_parm <- 0.1
nu_parm <- 1.5
lb <- rep(-Inf, n)
ub <- rep(1, n)
covmat <- fields::Matern(as.matrix(dist(locs)),
  range = range_parm,
  nu = nu_parm
)
L <- t(chol(covmat))
y <- as.vector(L %*% rnorm(n))
mask_cens <- y < ub
n_cens <- sum(mask_cens)
n_obs <- n - n_cens
rm(L)

# rearrange -------------------------
ind_cens <- which(mask_cens)
ind_obs <- setdiff(1:n, ind_cens)
new_order <- c(ind_obs, ind_cens)
y <- y[new_order]
locs <- locs[new_order, ]
covmat <- covmat[new_order, new_order]
lb <- lb[new_order]
ub <- ub[new_order]
mask_cens <- mask_cens[new_order]
ind_cens <- which(mask_cens)
ind_obs <- setdiff(1:n, ind_cens)

# three check plots
ind_test <- sample(which(mask_cens), 10, replace = FALSE)
check_obj <- ptmvn_check_converge(y, lb, ub, covmat,
  m_vec = seq(from = 10, to = 100, by = 10),
  ind_test = ind_test
)
first_mmt_ptmvn <- check_obj$pred
plot(first_mmt_ptmvn) + theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  axis.title.y = element_blank()
)
ggsave(
  paste0(
    "plots/mmt_cvg_ptmvn_lowdim.pdf"
  ),
  width = 6,
  height = 5
)

first_mmt_tmvn <- tmvn_check_converge(lb[mask_cens], ub[mask_cens],
  covmat[mask_cens, mask_cens],
  m_vec = seq(from = 10, to = 100, by = 10),
  ind_test = ind_test - n_obs
)
rownames(first_mmt_tmvn) <- paste("Loc", ind_test)
plot(first_mmt_tmvn) + theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  axis.title.y = element_blank()
)
ggsave(
  paste0(
    "plots/mmt_cvg_tmvn_lowdim.pdf"
  ),
  width = 6,
  height = 5
)

pred_err_tmvn <- first_mmt_tmvn - matrix(y[ind_test], nrow = length(ind_test),
                                         ncol = ncol(first_mmt_tmvn), byrow = F)
plot(pred_err_tmvn) + theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  axis.title.y = element_blank()
)  + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") 
ggsave(
  paste0(
    "plots/err_cvg_tmvn_lowdim.pdf"
  ),
  width = 6,
  height = 5
)
