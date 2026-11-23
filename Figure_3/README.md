# Difference-in-differences event-study figures

This directory contains panel data, Stata scripts and figures for four event-study analyses of antimicrobial use (AMU) and antimicrobial resistance (AMR) scores.

## Contents

```text
Figure_3/
├── README.md
├── code/       four analysis do-files
├── input/      four Stata datasets and their CSV companions
└── output/     four figures in PDF and PNG format, and Stata logs
```

| Analysis         | Script                       | Stata input                    | Output basename      | Event-time window                                           |
| ---------------- | ---------------------------- | ------------------------------ | -------------------- | ----------------------------------------------------------- |
| Animal AMR score | `code/Animal_AMR_score.do` | `input/Animal_AMR_score.dta` | `Animal_AMR_score` | −10 to +14                                                 |
| Animal AMU       | `code/Animal_AMU.do`       | `input/Animal_AMU.dta`       | `Animal_AMU`       | Full observed range                                         |
| Human AMR score  | `code/Human_AMR_score.do`  | `input/Human_AMR_score.dta`  | `Human_AMR_score`  | Full observed range; indicators are supplied in the dataset |
| Human AMU        | `code/Human_AMU.do`        | `input/Human_AMU.dta`        | `Human_AMU`        | −12 to +13; stacked DID bins the tails at −10 and +10     |

The scripts read the `.dta` files. Corresponding CSV files are provided for inspection; Stata metadata and numeric types are preserved in the `.dta` files.

## Run the analyses

The scripts target Stata 18 and require the user-written commands `ftools`, `reghdfe`, `did2s`, `event_plot` and `stackedev`, together with their dependencies. Install these commands using their respective package instructions before running. The scripts use the bundled `csdid` 1.82 and `drdid` 1.91 in `code/vendor/` ahead of installed versions. These files are from [friosavila/stpackages](https://github.com/friosavila/stpackages/tree/738defceb5413face7051754b259c82798c67e52) and include their MIT licence. Other package versions are not pinned.

In Stata, set the working directory to the `Figure_3` directory, then run:

```stata
do code/Animal_AMR_score.do
do code/Animal_AMU.do
do code/Human_AMR_score.do
do code/Human_AMU.do
```

 Each figure combines Callaway–Sant'Anna, Gardner DID-2S, stacked DID and TWFE estimates; event time −1 is the reference period.

All four scripts were run successfully using Stata/MP 18 to generate the included figures and logs. Event-study coefficient table rows in the new logs matched the original run's logged values for all four analyses.

## Contact

For questions on code, data, or reproducibility, contact the analysis team or the manuscript corresponding author.
