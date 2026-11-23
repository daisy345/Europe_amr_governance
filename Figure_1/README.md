# European antimicrobial use and resistance maps

This repository contains the inputs, R script and PDF outputs for four European choropleth maps of antimicrobial use (AMU) and antimicrobial resistance (AMR) scores in humans and animals.

## Repository contents

```text
Figure_1/
├── README.md
├── code/
│   └── Europe__hotmap.R
├── input/
│   ├── Human_AMR_score_hotmap.csv
│   ├── Human_AMU_hotmap.csv
│   ├── Animal_AMR_score_hotmap.csv
│   ├── Animal_AMU_hotmap.csv
│   └── geodata/
│       ├── Europe_Cleaned.shp
│       ├── Europe_Cleaned.shx
│       ├── Europe_Cleaned.dbf
│       └── Europe_Cleaned.prj
└── output/
    ├── Human_score.pdf
    ├── Human_AMU.pdf
    ├── Animal_AMR_score.pdf
    └── Animal_AMU.pdf
```

Keep the four boundary files together: the shapefile requires its accompanying `.shx`, `.dbf` and `.prj` files.

## Recreate the maps

Use R 4.6.2 or a compatible R installation with Cairo PDF support. Install the required packages once:

```r
install.packages(c("sf", "ggplot2", "scales"))
```

From the repository root, run:

```sh
Rscript code/Europe__hotmap.R
```

The script resolves input and output paths relative to its own location. It creates `output/` if needed and overwrites the four PDF files there. Each PDF is 10 × 8 inches.

| Input                           | Output                   | Colour scale             |
| ------------------------------- | ------------------------ | ------------------------ |
| `Human_AMR_score_hotmap.csv`  | `Human_score.pdf`      | 40–100; breaks every 20 |
| `Human_AMU_hotmap.csv`        | `Human_AMU.pdf`        | 0–50; breaks every 10   |
| `Animal_AMR_score_hotmap.csv` | `Animal_AMR_score.pdf` | 60–100; breaks every 20 |
| `Animal_AMU_hotmap.csv`       | `Animal_AMU.pdf`       | 0–200; breaks every 50  |

Missing values are shown in light grey.

## Data and reproducibility

Each CSV has one row per country, an `iso` country identifier and a numeric `value`. The script joins `iso` to `CNTR_ID` in the bundled boundary data; `country` and any other columns are not used in the join. It checks required columns and duplicate identifiers and warns when a CSV identifier is absent from the boundary data. The animal AMR score table contains precomputed values; the plotting script does not calculate the scores.

The maps were checked with R 4.6.2, sf 1.1.3, ggplot2 4.0.4 and scales 1.4.0. Different R, package or font versions may produce rendering differences.
