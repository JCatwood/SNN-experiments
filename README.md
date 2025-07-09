# Overview

This document provides guidance on installing the dependencies for as well as running the experiments that can replicate the results in **Scalable Sampling of Truncated Multivariate Normals Using Sequential Nearest-Neighbor Approximation**.

# Installing Dependencies

Most of our dependencies are on CRAN, with only one exception. The following R code can be used to installed the CRAN dependencies for running the experiments.
```
CRAN_pkg_names <- c("ggplot2", "fields", "scoringRules", "VeccTMVN", "sf", "spData", "GpGp", "nntmvn", "doParallel", "TruncatedNormal", "RColorBrewer", "autoimage", "mvtnorm", "RANN", "scales", "R.utils", "tidyr", "devtools")
for (pkg_name in CRAN_pkg_names) {
  install.packages(pkg_name, repos='http://cran.us.r-project.org')
}
```

`CensSpBayes` is the dependent R package that's not available on CRAN at this point. `CensSpBayes` is non-trivial to install. It is based on [`INLA`](https://www.r-inla.org/download-install), which requires several libraries that you may not be able to install through R. I used the following bash script on a freshly created Ubuntu VM to install the required libraries.
```
sudo apt-get install libudunits2-dev libssl-dev libgdal-dev libfontconfig-dev libgit2-dev libharfbuzz-dev libfribidi-dev libgsl-dev jags-dev
```

The following R code installs `CensSpBayes`. Notice that we made minor modification to the original implementation of [`CensSpBayes`](https://github.com/SumanM47/CensSpBayes) such that the MCMC samples can be returned.
```
library(devtools)
Git_pkg_name <- "CensSpBayes"
Git_pkg_link <- "JCatwood/CensSpBayes"
devtools::install_github(Git_pkg_link)
```

In fact, I could not install `CensSpBayes` on my Mac OS. The experiments involving `CensSpBayes` were run on a Linux system.

# Directories

## Section 3

Codes for producing figures/tables in Section 3 as well as the related figures/tables in the Appendix. All scripts can be run independently. 

- Figure 1 and Table 1, run `lowdim_SNN.R`, `lowdim_VMET.R`, and `lowdim_CSB.R`.

- Figure 2, run `time_compare.R` and `time_SNN.R`.

- Figure D.1, run `CSB_diagnosis.R`. Note that the successful installation of `CensSpBayes` R package is required.

## Section 4

Codes for producing figures/tables in Section 4 as well as the related figures/tables in the Appendix.

- Figure 3, run the first section of `performance_plot.R` with `score` set to `"RMSE"` and `order` set to `"maximin"`.
- Figure E.2, run the first section of `performance_plot.R` with `score` set to `"CRPS"` and `order` set to `"maximin"`.
- Figure E.3, run the first section of `performance_plot.R` with `score` set to `"RMSE"` and `order` set to `NULL`.
- Figure E.4, run the first section of `performance_plot.R` with `score` set to `"CRPS"` and `order` set to `NULL`.
- Figure E.5, run the second section of `performance_plot.R` with `score` set to `"RMSE"` and `"CRPS"`, sequentially.
- Figure 4, run `highdim_truth_known.R` with `use_snn_order` set to 2 for maximin ordering, `scene_ID` set to 1, and seed ID `k` set to 1. Run `highdim_truth_unknown.R` with `reorder` set to 1 for maximin ordering, `scene_ID` set to 1, and seed ID `k` set to 1.

The previous plots are based on the `mtd_cmp.csv` and the `m_cmp.csv` files. 

`mtd_cmp.csv` stores the (processed) outputs from running
- `highdim_truth_known.R` for 60 times, with `use_snn_order` set to 2, `scene_ID` sequentially set to 1, 2, and 3, and seed ID `k` sequentially set to $1, \ldots, 20$
- `highdim_truth_unknown.R` for 60 times, with `reorder` set to 1, `scene_ID` sequentially set to 1, 2, and 3, and seed ID `k` sequentially set to $1, \ldots, 20$
- `highdim_truth_known.R` for 60 times, with `use_snn_order` set to 0, `scene_ID` sequentially set to 1, 2, and 3, and seed ID `k` sequentially set to $1, \ldots, 20$
- `highdim_truth_unknown.R` for 60 times, with `reorder` set to 0, `scene_ID` sequentially set to 1, 2, and 3, and seed ID `k` sequentially set to $1, \ldots, 20$

`m_cmp.csv` stores the (processed) outputs from running
- `m_compare_truth_known.R` for 180 times, with `use_snn_order` sequentially set to 0, 1, and 2; `scene_ID` sequentially set to 1, 2, and 3; and seed ID `k` sequentially set to $1, \ldots, 20$

`performance_table.R` is legacy, which is not used for any presentation in the current version.

## Section 5

Figures 5 and 6 in Section 5 can be reproduced by running `PCE_modeling.R`, where the control variables, `run_parm_est`, `run_VeccTMVN`, `run_TN`, `run_CB`, `run_seq_Vecc`, `run_seq_Vecc_all`, and `run_CB_all` should be set to `TRUE`.
 
Certain components of this script can take relatively long time

- `run_parm_est` is the switch for covariance parameter estimation, which took around 10 hours. Its result was store in `results/PCE_modeling.RData`. You can reuse the previous covariance estimation by setting `run_parm_est` to `FALSE` 
- `run_CB` and `run_CB_all` are the switches for using the CSB method for the inference over Texas and the entire US, respectively, the latter of which may take 3 hours
- All other methods are relatively fast












