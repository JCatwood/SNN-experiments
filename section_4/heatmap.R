rm(list = ls())

library(fields)

m <- 30
k <- 1
scene_ID <- 1
use_snn_order <- 0
source("../utils/data_simulation.R")

y_true <- read.table(paste0("plots/sim_data_scene", scene_ID, "_true.csv"),
  header = F, row.names = NULL
)[, 1]
y_SNN <- read.table(
  paste0(
    "plots/samp_cmp_known_SNN_scene", scene_ID, "_m", m,
    "_order", use_snn_order, ".csv"
  ),
  header = F, row.names = NULL
)[, 1]
y_VT <- read.table(
  paste0(
    "plots/samp_cmp_known_VT_scene",
    scene_ID, ".csv"
  ),
  header = F, row.names = NULL
)[, 1]
y_TN <- read.table(
  paste0(
    "plots/samp_cmp_known_TN_scene",
    scene_ID, ".csv"
  ),
  header = F, row.names = NULL
)[, 1]
y_CB <- read.table(
  paste0(
    "plots/samp_cmp_known_CB_scene",
    scene_ID, ".csv"
  ),
  header = F, row.names = NULL
)[, 1]

zlim <- range(y_true, y_SNN, y_TN, y_VT, y_CB)

pdf(
  file = paste0("plots/sim_data_scene", scene_ID, "_true.pdf"),
  width = 5, height = 5
)
fields::image.plot(matrix(y_true, sqrt(n + n_test), sqrt(n + n_test)), zlim = zlim)
dev.off()

pdf(
  file = paste0(
    "plots/samp_cmp_known_SNN_scene", scene_ID, "_m", m,
    "_order", use_snn_order, ".pdf"
  ),
  width = 5, height = 5
)
fields::image.plot(matrix(y_SNN, sqrt(n + n_test), sqrt(n + n_test)), zlim = zlim)
dev.off()

pdf(
  file = paste0("plots/samp_cmp_known_VT_scene", scene_ID, ".pdf"),
  width = 5, height = 5
)
fields::image.plot(matrix(y_VT, sqrt(n + n_test), sqrt(n + n_test)), zlim = zlim)
dev.off()

pdf(
  file = paste0("plots/samp_cmp_known_TN_scene", scene_ID, ".pdf"),
  width = 5, height = 5
)
fields::image.plot(matrix(y_TN, sqrt(n + n_test), sqrt(n + n_test)), zlim = zlim)
dev.off()

pdf(
  file = paste0("plots/samp_cmp_known_CB_scene", scene_ID, ".pdf"),
  width = 5, height = 5
)
fields::image.plot(matrix(y_CB, sqrt(n + n_test), sqrt(n + n_test)), zlim = zlim)
dev.off()
