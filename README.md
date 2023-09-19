# SpatialProbabilisticMetaweb


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8350065.svg)](https://doi.org/10.5281/zenodo.8350065)

Code repo for the manuscript *Spatially explicit predictions of food web structure from regional level data* available as a preprint on [EcoEvoRxiv](https://ecoevorxiv.org/repository/view/5941/) (manuscript preparation on [GitHub](https://github.com/PoisotLab/ms_spatial_metaweb)).

## Scripts

All scripts are in the top level folder. [`main.jl`](main.jl) contains instructions to run all analysis scripts. Briefly:

- `01_...` and other numbered scripts contain the main analysis steps.
- `A0_...` and similar scripts define useful functions for the analyses.
- `x1_...` and further scripts are preparation scripts that produce elements available on the repo (e.g. in the data/input/ folder).

## Data

Most data files are too large to be hosted here, especially intermediate analysis files, so *do not expect to be able to run the scripts directly yet*. Almost all files in `data` can be downloaded or produced using the scripts though. The long-term goal is to have everything necessary to reproduce a smaller-scale example of the analyses.

## Other folders

- `figures`: contains all the figures produced in the analysis steps
- `jobs`: contains the job scripts to run the resource-intensive scripts on compute clusters
- `xtras`: contains previous scripts (version-controlled) and the files to run a smaller-scale analysis (not version-controlled yet)
