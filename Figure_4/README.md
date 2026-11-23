# Reproduction of Figure 4

This repository contains processed inputs and an R script for two figures of European mean human antimicrobial resistance (HAMR) scores across SSP pathways and policy settings.

```text
Figure_4/
├── README.md
├── input/
│   ├── historical_analysis_data_2007_2021.csv
│   └── future_scenario_covariates_2022_2050.csv
├── code/
│   └── Figure_4.R
└── output/
    ├── Figure_4a.pdf
    └── Figure_4b.pdf
```

`Figure_4a.pdf` shows historical and projected trajectories through 2050, including 95% confidence intervals for projections. `Figure_4b.pdf` shows the change between 2021 and 2050. The inputs are processed data; upstream raw datasets are not included.

## Reproduce the figures

Use R 4.6.2 or a compatible installation with Cairo PDF support, the recommended package `nlme`, and `ggplot2`. Install `ggplot2` if necessary:

```r
install.packages("ggplot2")
```

From the repository root, run:

```sh
Rscript code/Figure_4.R
```

The script fits the model and produces both figures using a temporary directory for intermediate files. It writes `Figure_4a.pdf` and `Figure_4b.pdf` to `output/`, overwriting files with those names. The workflow was tested on Windows with R 4.6.2, nlme 3.1-171 and ggplot2 4.0.4; results and rendering may vary with the software environment.

## Analysis

The script validates a country random-intercept model on 2019–2021 after fitting it to 2007–2018 data, then fits the final model to 2007–2021. It projects outcomes under SSP1, SSP2, SSP3 and SSP5 and externally specified policy-effect schedules. The future confidence intervals include fixed-effect and posterior country random-intercept uncertainty, but not uncertainty in future covariates or the policy-effect schedule.

Before public release, document the upstream data sources, variable definitions, units and citations, confirm redistribution rights for the processed data, and choose appropriate data and code licences.

## Contact

For questions on code, data, or reproducibility, contact the analysis team or the manuscript corresponding author.
