# Summary

Experiments for **Scalable Sampling of Spatial Truncated Multivariate Normals Using Sequential Nearest-Neighbor Approximation**.

## Dependencies

```r
install.packages(c("nntmvn", "TruncatedNormal", "VeccTMVN", "GpGp", "fields",
  "scoringRules", "RANN", "ggplot2", "dplyr", "tidyr", "lhs",
  "doParallel", "R.utils", "sf", "spData", "mvtnorm", "knitr"))
```

`CensSpBayes` is required to reproduce all numerical results in the paper but optional for experiments not involving the CSB method. Its modified version at `JCatwood/CensSpBayes` returns MCMC samples.

## Section 3

Run from `section_3`:
```bash
mkdir -p plots
mkdir -p results
```
For Figure 1 and Table 1, run
```
Rscript lowdim_SNN.R
Rscript lowdim_VMET.R
Rscript lowdim_CSB.R
```
Notice that `lowdim_CSB.R` requires the existence of the `CensSpBayes` package. For Figure 2, run
```bash
Rscript converge_analysis.R 1000
```
Outputs are `results/local_moments_lowdim.csv`, a reproducibility RDS file, and `plots/local_moments_lowdim.pdf`. Each neighborhood is sampled with MET. Means and unbiased sample variances use the same draws; Monte Carlo standard errors are saved. Ten reproducible target indices are used. 

For Figure 3, run
```bash
Rscript time_compare.R
Rscript time_SNN.R
```

## Section 4

Run from `section_4`:
```bash
mkdir -p results/logs
mkdir -p plots
```
Run the following script for data simulation once (not in parallel):
```bash
for scene in 1 2 3; do
  Rscript -e "scene_ID <- $scene; source('../utils/data_simulation.R')"
done
```
Delete `section_4/data/` whenever the scenario definitions change. For the following experiment codes, users should first set the `seed` variable (chosen from 1 to 20). For example
```bash
export seed=1
```
Then run
```bash
for scene in 1 2 3; do
  Rscript SNN_m_cmp_known.R "$seed" "$scene" 2 > "results/logs/snn_${scene}_${seed}.log" 2>&1
  Rscript VMET_m_cmp_known.R "$seed" "$scene" > "results/logs/vmet_${scene}_${seed}.log" 2>&1
  Rscript MET_known.R "$seed" "$scene" > "results/logs/met_${scene}_${seed}.log" 2>&1
done
```
for the known covariance simulations. Here, the `2` in `SNN_m_cmp_known.R "$seed" "$scene" 2` indicates maximin ordering. Run from `section_4`:
```bash
for scene in 1; do
  Rscript CSB_unknown.R "$seed" > "results/logs/CSB_unknown_${scene}_${seed}.log" 2>&1
  Rscript SNN_unknown.R "$seed" > "results/logs/SNN_unknown_${scene}_${seed}.log" 2>&1
done
```
for the unknown covariance simulations. Notice that `CSB_unknown.R` requires the existence of the `CensSpBayes` package. Run from `section_4`:
```bash
for scene in 1; do
  Rscript SNN_m_cmp_known.R "$seed" "$scene" 0 > "results/logs/snn_coord_${scene}_${seed}.log" 2>&1
  Rscript SNN_m_cmp_known.R "$seed" "$scene" 1 > "results/logs/snn_vardesc_${scene}_${seed}.log" 2>&1
done
```
to get the results for two other ordering schemes considered for the SNN method. To produce Figures 4 and 8, run `performance_plot.R`. For Table 4, run `performance_table.R`. Notice that 600 seconds were added to `MET` and `VMET` because their first partition is needed and partitions are applied when a timeout is encountered. To produce Figure 5, run `sample_heatmap.R`, which uses a saved sample from CSB and newly generates samples for SNN, MET, VMET. This script does not require the `CensSpBayes` package.

## Section 5

For Figures 6, run from `section_5`:

```bash
for seed in {1..10}; do
  Rscript PCE_validate.R TX "$seed" 50 SNN,VMET,MET,LOD
  Rscript PCE_validate.R US "$seed" 50 SNN,LOD
done
Rscript PCE_figures.R
```

Add `CSB` to either method list on a machine with that package. Regional MET/VMET may still be expensive. The national default omits them because a full training TMVN may be infeasible.

The split holds out about 20% of distinct longitude/latitude pairs, so repeated visits to the same well stay together. All response standardization and covariance fitting use training records only. A covariance fit is recomputed per split.

Gaussian methods share a local marginal prediction rule with 50 nearest training responses in ARD-scaled coordinates and include residual Gaussian uncertainty. CSB uses its own spatial-only predictive model. Outcomes are scored with Brier score for predicted censoring probability where detection limits are known.

Results are saved in `results/`.

For Figure 7, in `PCE_modeling.R`, enable `run_seq_Vecc_all`. Leave `run_parm_est=FALSE` to use the supplied full-data fit for these descriptive maps. 
```bash
Rscript PCE_modeling.R
Rscript PCE_figures.R
```
The national plots visualize predicted posterior means and their standard errors at the censored records.
