# Antimicrobial use, resistance and AMR-associated death trends

This directory contains processed inputs, R plotting scripts and figures for trends in human and animal antimicrobial use (AMU), human and animal antimicrobial resistance (AMR) scores, and AMR-associated deaths.

## Contents

```text
Figure_2/
├── README.md
├── code/       five R plotting scripts
├── input/      five Joinpoint exports and accompanying CSV or Excel tables
└── output/     five figures in PDF and PNG format
```

| Analysis | Script in `code/` | Joinpoint export in `input/` | Outputs in `output/` |
| --- | --- | --- | --- |
| Human AMU | `human_AMU.R` | `human AMU.txt` | `Human_AMU.pdf`, `Human_AMU.png` |
| Animal AMU | `animal_AMU.R` | `animal AMU.txt` | `Animal_AMU.pdf`, `Animal_AMU.png` |
| Human AMR score | `human_AMR_score.R` | `human AMR score.txt` | `Human_AMR_score.pdf`, `human_AMR_score.png` |
| Animal AMR score | `animal_AMR_score.R` | `animal AMR score.txt` | `Animal_AMR_score.pdf`, `Animal_AMR_score.png` |
| AMR-associated deaths | `AMR-associated death.R` | `AMR-associated death.txt` | `AMR-associated_death.pdf`, `AMR-associated_death.png` |

## Analysis

Joinpoint regression was performed using Joinpoint Regression Software (version 5.4.0; National Cancer Institute, USA) to identify changes in temporal trends in human and animal antimicrobial resistance (AMR) scores and antimicrobial use (AMU).

The statistical analyses were conducted in Joinpoint Regression Software. The resulting trends were subsequently replotted using R code to produce the final figures.

The supplied text files contain observed values, fitted values (`Model`), joinpoint flags and annual percent changes (APC). The scripts plot these exported results; they do not fit Joinpoint models or calculate APCs or statistical significance. APC labels and some joinpoint annotations are fixed in the scripts for these results and should be checked if the inputs are changed.

The CSV files for AMU and Excel files for AMR scores and deaths are accompanying tables; the R scripts do not read them. The original Joinpoint project, full model settings and complete upstream data-processing workflow are not included. The plots can be recreated from the supplied exports, but the underlying statistical analyses cannot be independently rerun from this directory alone.

## Reproduce the figures

Use R 4.6.2 or a compatible R installation. The scripts use base R and require no additional packages. From the `Figure_2` directory, run:

```sh
Rscript "code/human_AMU.R"
Rscript "code/animal_AMU.R"
Rscript "code/human_AMR_score.R"
Rscript "code/animal_AMR_score.R"
Rscript "code/AMR-associated death.R"
```

Each script resolves paths relative to its own location, creates `output/` if needed, and overwrites the corresponding PDF and PNG. Each figure is 7.1 × 5 inches and each PNG is 600 dpi. All five scripts were run successfully on Windows with R 4.6.2; rendering and fonts may vary with the software environment.

Before public release, document the upstream data sources, variable definitions, units and citations, confirm redistribution rights for the processed data, and choose appropriate data and code licences.

## Contact

For questions on code, data, or reproducibility, contact the analysis team or the manuscript corresponding author.
