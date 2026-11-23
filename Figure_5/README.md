# Reproduction of Figure 5

This repository contains processed inputs and an R script for mortality projections under four SSP pathways and four human antimicrobial resistance (HAMR) scenarios.

```text
Figure_5/
├── README.md
├── input/
│   ├── historical_panel.csv
│   ├── future_covariates.csv
│   └── hamr_scenario_trajectories.csv
├── code/
│   └── Figure_5.R
└── output/
    ├── covariate_forecast_backtest.pdf
    ├── death_change_2050_vs_2019.pdf
    └── death_predictions_2007_2050i.pdf
```

The figures show the mortality-model backtest, projected percentage changes in deaths from 2019 to 2050, and predicted death trajectories. The inputs contain the historical country-year panel, future covariates and country-level HAMR scenarios; upstream raw data and the model used to create the HAMR scenarios are not included.

## Reproduce the figures

Use R 4.6.2 or a compatible R installation with Cairo graphics support. Install the required packages if needed:

```r
install.packages(c("data.table", "fixest", "ggplot2", "scales"))
```

From the repository root, run:

```sh
Rscript code/Figure_5.R
```

The script fits and validates the mortality model, predicts deaths and generates the three figures. Intermediate files are created in a temporary directory and removed when the script finishes. The three existing files in `output/` are overwritten. The workflow was tested with R 4.6.2, data.table 1.18.6.1, fixest 0.14.2, ggplot2 4.0.4 and scales 1.4.0 on Windows. The figures use the Arial font; rendering may differ in other environments.

## Analysis

The backtest trains on data through 2018 and evaluates predictions for 2019–2021 across HAMR lag orders 0–6. It selects the adjusted model with the lowest log-scale test RMSE and refits it on the available historical data. With the supplied inputs, the selected maximum lag is four years and test RMSE is approximately 0.06946. Predictions begin in 2011 because four earlier years of HAMR observations are needed.

The trajectory bands sum country-level confidence limits; they do not include uncertainty in the upstream covariates or HAMR scenarios. The 2050 change figure uses the corresponding SSP's predicted 2019 HAMR-current total as its baseline.

Before public release, document the upstream data sources, variable definitions, units and citations, confirm redistribution rights for the processed inputs, and select appropriate data and code licences.

## Contact

For questions on code, data, or reproducibility, contact the analysis team or the manuscript corresponding author.
