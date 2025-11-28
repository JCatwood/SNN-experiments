library(GpGp)
library(VeccTMVN)
library(TruncatedNormal)
library(ggplot2)
library(nntmvn)
# highdim ------------------------------
## simulation settings ------------------------------
rm(list = ls())
set.seed(123)
m <- 30 # number of nearest neighbors
n_samp <- 50 # samples generated for posterior inference
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  k <- as.integer(args[1]) # k is the index for GP realizations
  scene_ID <- as.integer(args[2])
} else {
  k <- 1
  scene_ID <- 3
}

## data simulation ----------------------
source("../utils/data_simulation.R")
y <- y_list[[k]]
mask_cens <- (y < cens_ub) & (y > cens_lb)

## convergence check ------------------------
check_obj <- ptmvn_check_converge(y, cens_lb , cens_ub, covmat,
                                  m_vec = seq(from = 10, to = 100, by = 10))
pred_err <- check_obj$error
plot(pred_err) + theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  axis.title.y = element_blank()
)
ggsave(
  paste0(
    "plots/converg_analysis_scene_", scene_ID, ".pdf"
  ),
  width = 6,
  height = 5
)