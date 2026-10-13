# RI-CLPM analyses of desistance correlates: R code (Pathways to Desistance, ICPSR 29961)

These R scripts build the annual person-year analysis file from the Pathways to Desistance public-use data and reproduce every table, figure and simulation reported in the manuscript and its Online Resource 1. The main models are random-intercept cross-lagged panel models (RI-CLPMs) fitted in lavaan for ages 16 to 23, with maximum-likelihood dynamic panel models (dpm), multivariate RI-CLPMs and three simulation studies as checks.

## Data (not included)

The data are not part of this repository. The ICPSR terms of use do not allow redistribution. Registered ICPSR users can download them here:

Mulvey, E. P. (2016). *Research on Pathways to Desistance [Maricopa County, AZ and Philadelphia County, PA]: Subject measures, 2000–2010* (ICPSR 29961; Version 2). Inter-university Consortium for Political and Social Research. https://doi.org/10.3886/ICPSR29961.v2

Download the R version of the data (`29961-0001-Data.rda` to `29961-0011-Data.rda`). The `reference_outputs/` folder in this repository holds only summary estimates (coefficients, intervals, fit indices, simulation results), never respondent-level data.

## Folder layout

Put this repository inside the folder that holds the ICPSR download:

```
ICPSR_29961/
  DS0001/29961-0001-Data.rda
  ...
  DS0011/29961-0011-Data.rda
  pathways-riclpm-replication/      <- this repository
    00_RUN_ALL.R
    R/
    reference_outputs/
```

## How to run

1. Install the packages once:
   ```r
   install.packages(c("lavaan", "dpm", "panelr", "mice", "sandwich",
                      "lmtest", "fixest", "ggplot2"))
   ```
2. Open `00_RUN_ALL.R` in RStudio and click Source (or run `source("pathways-riclpm-replication/00_RUN_ALL.R")` from the `ICPSR_29961` folder). The script finds its own folder; if it cannot, edit `pkg_dir` near the top.
3. The full run took 63 minutes on a Windows laptop (R 4.6.1). Most of the time goes to `05_riclpm_main.R` and `09_multivariate_riclpm.R`. lavaan warns that some pairs of ages have 0% coverage. That is expected, because no youth was interviewed at every age from 16 to 23.

Everything is written inside `ICPSR_29961`:

| Output | Written by |
|---|---|
| `panel_desistance.rds`, `panel_6month.rds` | 01 |
| `output_checks/` | 03, 04 |
| `output_main/` | 05, 06, 09 |
| `figures/` (Fig1 and Fig2 as 600 dpi TIFF and EPS) | 07 |
| `verification_report.txt` | 08 |

If an older `panel_desistance.rds` exists, it is copied to `panel_desistance_before_2026-10-12.rds` before 01 replaces it.

## Scripts

| Script | What it does |
|---|---|
| `00_RUN_ALL.R` | Runs the scripts below in order and lists any step that failed. |
| `R/01_data_prep.R` | Builds the annual person-year panel (ages from baseline age plus years since baseline). |
| `R/03_robustness_checks.R` | Earlier estimators: lagged hybrid models with 20 multiple imputations, fixed-effects Poisson, a positive control, within-person variation and first RI-CLPMs. |
| `R/04_simulation_lag_artifacts.R` | Simulations 1 and 2: lagged estimates when the only link is within the same period, and calibrated null cross-lags. |
| `R/05_riclpm_main.R` | The 18 RI-CLPMs (14 scales and 4 social bonds) and their checks: age bands against one cross-lag for all ages, time in secure facilities, same-age directed paths, childhood factors, dynamic panel models and random slopes. |
| `R/06_simulation_common_state.R` | Simulation 3: a shared, persistent state with no causal effect in either direction. |
| `R/09_multivariate_riclpm.R` | Four multivariate RI-CLPMs (cognitive, maturation, bonds and a cross-domain set), so each cross-lag holds the related measures constant. |
| `R/07_figures.R` | Figures 1 and 2. |
| `R/08_verify_against_reference.R` | Compares a new run with `reference_outputs/` (MATCH, CLOSE or DIFFERS for each table) and recomputes 85 numbers quoted in the manuscript (PASS or FAIL). |

The numbering skips 02. That script belonged to an earlier version of the analysis and is not needed.

## Reproducibility

Every script sets a fixed random seed. The tables in `reference_outputs/` come from a run with R 4.6.1, lavaan 0.7-2, dpm 1.3.0, panelr 1.0.1, mice 3.19.0, fixest 0.14.2, sandwich 3.1-3, lmtest 0.9-40 and ggplot2 4.0.3. Against those tables, that run gave 34 of 34 tables MATCH and 85 of 85 checks PASS.

An earlier run with R 4.3.3 and lavaan 0.6-17 gave identical estimates. Two things differed across versions: the robust fit indices of the four social-bond models, and the simulation results (the random draws differ). Expect the same with other versions.

## Citation

Please cite the archived version of this code (see `CITATION.cff` or the Zenodo record) and the ICPSR data above.

## License

The code is released under the MIT License (see `LICENSE`).
