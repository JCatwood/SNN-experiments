# SNN manuscript revision

Experiments for **Scalable Sampling of Spatial Truncated Multivariate Normals Using Sequential Nearest-Neighbor Approximation**.

## Dependencies

```r
install.packages(c("nntmvn", "TruncatedNormal", "VeccTMVN", "GpGp", "fields",
  "scoringRules", "RANN", "ggplot2", "dplyr", "tidyr", "scales", "lhs",
  "doParallel", "R.utils", "sf", "spData", "RColorBrewer", "autoimage"))
```

`CensSpBayes` is required to reproduce all numerical results in the paper but optional for experiments not involving the CSB method. Its modified version at `JCatwood/CensSpBayes` returns MCMC samples.

## Local mean and variance diagnostics

Run from `section_3`:

```bash
Rscript converge_analysis.R 1000
# Optional full-neighborhood MET reference (can be expensive):
Rscript converge_analysis.R 1000 true
```

Outputs are `results/local_moments_lowdim.csv`, a reproducibility RDS file, and `plots/local_moments_lowdim.pdf`. Each neighborhood is sampled with MET. Means and unbiased sample variances use the same draws; Monte Carlo standard errors are saved. Ten reproducible target indices are used. 
If `Rscript converge_analysis.R 1000 true` was run, full-reference rows are saved in the CSV/RDS and intentionally omitted from the m=10,...,100 plot; inspect them to assess reference error.

## Simulation study

The simulation generates 10,000 responses and censors those below 1. Scenario 2 uses a normalized random Gram covariance with a 0.1 nugget. Scenario 3 has four-dimensional Latin-hypercube inputs.

For the following codes, users should first set the `seed` variable (chosen from 1 to 20). Run from `section_4`:
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
To produce Figures 4 and 8, run `performance_plot.R`. For Table 4, run `performance_table.R`.

To produce Figure 5, run `sample_heatmap.R`, which uses a saved sample from CSB and newly generates samples for SNN, MET, VMET. This script does not require the `CensSpBayes` package.

## PCE validation and Figure 6

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

## National summaries and Figure 7

In `PCE_modeling.R`, enable `run_seq_Vecc_all`. Leave `run_parm_est=FALSE` to use the supplied full-data fit for these descriptive maps. 

```bash
Rscript PCE_modeling.R
Rscript PCE_figures.R
```

The national plots visualize predicted posterior means and their standard errors at the censored records.
