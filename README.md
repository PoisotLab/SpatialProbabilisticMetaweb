# SpatialProbabilisticMetaweb


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8350065.svg)](https://doi.org/10.5281/zenodo.8350065)

Code repo for the manuscript *Spatially explicit predictions of food web structure from regional level data* available as a preprint on [EcoEvoRxiv](https://ecoevorxiv.org/repository/view/5941/) (manuscript preparation on [GitHub](https://github.com/PoisotLab/ms_spatial_metaweb)).

## Scripts

All scripts for the main analyses are numbered and in the top level folder. [`00_main.jl`](00_main.jl) contains instructions to run all analysis steps. Additional scripts are in the `scripts/` folder, where:

- `lib/`  contains useful functions for the analyses.
- `prep/`  contains preparation scripts that produce elements available on the repo (e.g. in the `data/input/` folder) and required to run the main analyses.

## Data

Most data files are too large to be hosted here, especially intermediate analysis files, so *do not expect to be able to run the top-level scripts directly*. Almost all files in `data/` can be downloaded or produced using the scripts in `scripts/prep/`. Note that running the analyses at full-scale is resource-intensive and required the use of high-memory compute clusters.

## Other folders

- `figures`: contains all the figures produced in the analysis steps
- `jobs`: contains the job scripts to run the resource-intensive scripts on compute clusters
- `xtras`: contains previous scripts (version-controlled) and the files to run a smaller-scale analysis (not version-controlled yet)
