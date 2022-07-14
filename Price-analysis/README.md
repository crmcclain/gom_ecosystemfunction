Deep-sea BEF price analysis
================

This project contains the script to perform all Price analyses on
woodfall communities.

`01_price-data-load.R` performs:

1)  data import and cleaning from raw-data files ‘logdiversity.rds’

2)  establishing site order and performing pairwise comparisons

3)  performing data transformations on price components to relative
    importance, relative contribution, and normalized scales for
    analysis and summary

`02_price_work.R` subsequently:

1)  outputs summaries

2)  creates relative importance plot

3)  creates relative contribution plot

4)  loads generalized additive models (GAMs) (all code to rerun is found
    in `run-GAM-script.R`)

5)  creates plot of predicted effects of deltaN on Price components

6)  combines all plots for publication output into ./figures/ folder

The directory structure is:

    .
    ├── Price-analysis
    │   ├── code
    │   │   ├── 01_price-data-load.R
    │   │   ├── 02_price_work.R
    │   │   └── run-GAM-script.R
    │   ├── data
    │   │   ├── derived-data
    │   │   │   ├── models
    │   │   │   │   ├── CDE_null.rds
    │   │   │   │   ├── CDE_n_brm.rds
    │   │   │   │   ├── CDE_SIE.rds
    │   │   │   │   ├── CDE_SRE.rds
    │   │   │   │   ├── SIE.G_null.rds
    │   │   │   │   ├── SIE.G_n_brm.rds
    │   │   │   │   ├── SIE.L_null.rds
    │   │   │   │   ├── SIE.L_n_brm.rds
    │   │   │   │   ├── SRE.G_null.rds
    │   │   │   │   ├── SRE.G_n_brm.rds
    │   │   │   │   ├── SRE.L_null.rds
    │   │   │   │   └── SRE.L_n_brm.rds
    │   │   │   ├── N1000Df.rds
    │   │   │   ├── N500Df.rds
    │   │   │   ├── pricePartitions.rds
    │   │   │   └── pricePartitionsN.rds
    │   │   └── raw-data
    │   │       ├── logdiversity.rds
    │   │       └── logdiversity_drop.rds
    │   ├── figures
    │   │   └── Price_nEffects.pdf
    │   ├── README.md
    │   └── README.Rmd
