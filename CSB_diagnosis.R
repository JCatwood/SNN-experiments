#####################################################################
#
# Produces the MCMC diagnosis plots for CensSpBayes
#
#####################################################################

library(ggplot2)

rm(list = ls())
set.seed(123)

load("./results/sim_data_CB_scene2_rep1.RData")
ind <- sample(1:nrow(y_samp_CB$Y.pred.samp), 6, replace = FALSE)
j <- 1
for (i in ind) {
  mydf <- data.frame(
    x = 1:ncol(y_samp_CB$Y.pred.samp),
    y = y_samp_CB$Y.pred.samp[i, ]
  )
  ggplot(data = mydf, mapping = aes(x = x, y = y)) +
    geom_point(size = 1) +
    scale_x_continuous(name = "MCMC sample index after burn-in and thinning") +
    scale_y_continuous(name = paste0(i, "-th censored response")) +
    theme(
      text = element_text(size = 14),
      plot.title = element_blank()
    )
  if (!file.exists("plots")) {
    dir.create("plots")
  }
  ggsave(
    paste0("plots/CSB_diagnosis_", j, ".pdf"),
    width = 5,
    height = 5
  )
  j <- j + 1
}
