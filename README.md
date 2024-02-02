# SpatialProbabilisticMetaweb


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8350065.svg)](https://doi.org/10.5281/zenodo.8350065)

Code repo for the manuscript *Spatially explicit predictions of food web structure from regional level data* available as a preprint on [EcoEvoRxiv](https://ecoevorxiv.org/repository/view/5941/) (manuscript preparation on [GitHub](https://github.com/PoisotLab/ms_spatial_metaweb)). This repository provides code for both the full-scale analyses and a reproducible minimal example (see instructions below).

## Scripts

All scripts for the main analyses are numbered and in the top level folder. [`00_main.jl`](00_main.jl) contains instructions to run all analysis steps. Additional scripts are in the `scripts/` folder, where:

- `lib/`  contains useful functions for the analyses.
- `prep/`  contains preparation scripts that produce elements available on the repo (e.g. in the `data/input/` folder) and required to run the main analyses.

## Folder organization

- `data`: contains the input data and result layers. Intermediate outputs of full-scale analyses are also exported here but are not version-controlled.
- `figures`: contains all the figures produced in the analysis steps.
- `jobs`: contains the job scripts to run the resource-intensive scripts on compute clusters.
- `xtras`: contains the input files to run a reproducible minimal example. Results (version-controlled) and intermediate outputs are also exported here.

## Reproducibility

### Install Julia

We recommend installing Julia using the cross-platform installer [juliaup](https://github.com/JuliaLang/juliaup), which makes it easy to install a specific Julia version. We used Julia 1.9.1 for this project.

1. Follow the [juliaup installation instructions](https://github.com/JuliaLang/juliaup).
2. Install Julia 1.9.1 through juliaup in a terminal:
```
juliaup add 1.9.1
```

### Setting up the environment

We used the Julia package manager to track the packages and versions used in the [Project.toml](Project.toml) and [Manifest.toml](Manifest.toml) files. To set up your environment similarly:

1. Clone this repository and launch Julia from a terminal in the top level folder. Make sure to always launch Julia this way using the project.

```
julia +1.9.1 --project
```

2. Install packages 

```julia
using Pkg; Pkg.instantiate()
```

### Running the scripts

- *All analyses are reproducible using the scripts in this repo*. However, not all of them are directly reproducible from the repository because of file size restrictions and required computations (see the following table). By direct reproducibility we mean cloning the repo, setting up Julia as detailed above, then running the script *alone and as-is* (without re-running previous scripts).
- **Minimal example**: Scripts are directly reproducible with the minimal example. Leave the scripts as-is or set `CAN = false` before running. You can directly run scripts `06`, `07`, etc.
- **Full-scale analyses**: Plotting scripts are directly reproducible at full scale. All other analysis steps are scripted but as a whole they require high resources (most steps were run on high memory clusters). To run analyses at full scale, manually set `CAN = true` before running the script.
- Many steps use [threads](https://docs.julialang.org/en/v1/manual/multi-threading/) (if available) to run analyses faster. We recommend doing so if you have multiple cores available, but note that it will require more memory. For example, to start Julia with 8 threads available: `julia +1.9.1 --project --threads 8`
- The following table shows the reproducibility of all scripts. Timing for the minimal example was evaluated on a laptop with 8 GB of RAM and 8 cores (using 8 threads). Running [00_main.jl](00_main.jl), which re-runs all other scripts, took ~ 1h 20 minutes.

| Script | Directly reproducible for minimal example (10 arc-min) | Directly reproducible for CAN (2.5 arc-min) |
| ---- | ---- | ---- |
|  | Leave as-is or set `CAN = false` | Set `CAN = true` |
| **Main** |  |  |
| [00_main.jl](00_main.jl) | ✅ Yes, re-runs everything in ~ 1h 20 min | 🚫 Not recommended |
|  |  |  |
| **Data preparation** |  |  |
| [01_get_occurrences.jl](01_get_occurrences.jl) | ✅ Yes | ✅ Same as minimal |
| [02_get_absences.jl](02_get_absences.jl) | ✅ Yes<br>⏰ ~ 25 minutes | ✔️ Data available<br>:warning: Requires several hours |
| [03_generate_sdms.jl](03_generate_sdms.jl) | ✅ Yes<br>⏰ ~ 12 minutes | ⤴ Requires previous<br>:warning: Memory intensive |
|  |  |  |
| **Assembling results** |  |  |
| [04_aggregate_sdms.jl](04_aggregate_sdms.jl) | ✅ Yes<br>⤴ Will re-run previous script once | ⤴ Requires previous<br>:warning: Memory intensive |
| [05_assemble_networks.jl](05_assemble_networks.jl) | ✅ Yes | 🚫 Requires high memory and long computations |
|  |  |  |
| **Analysis** |  |  |
| [06_get_species_lcbd.jl](06_get_species_lcbd.jl) | ✅ Yes<br> | ⤴ Requires 04<br>:warning: Memory intensive |
| [07_get_network_lcbd.jl](07_get_network_lcbd.jl) | ✅ Yes<br> | 🚫 Requires high memory and long computations |
| [09_get_network_measures.jl](09_get_network_measures.jl) | ✅ Yes<br> | 🚫 Requires high memory and long computations |
| [11_get_ecoregions_measures.jl](11_get_ecoregions_measures.jl) | ✅ Yes | ✅  Yes |
|  |  |  |
| **Plotting results** |  |  |
| [08_plot_lcbd_results.jl](08_plot_lcbd_results.jl) | ✅ Yes | ✅ Yes |
| [10_plot_network_measures.jl](10_plot_network_measures.jl) | ✅ Yes | ✅ Yes |
| [12_plot_ecoregions.jl](12_plot_ecoregions.jl) | ✅ Yes | ✅ Yes |
|  |  |  |
| **Motifs analysis** |  |  |
| ⚠ This analysis is especially intensive ⚠ |  |  |
| [13_get_motifs.jl](13_get_motifs.jl) | ✅ Yes<br>:warning: ~ 25 mins per motif, only for New-Brunswick | 🚫🚫🚫 Requires job arrays |
| [14_assemble_motifs.jl](14_assemble_motifs.jl) | ✅ Yes | 🚫 Requires previous |
| [15_plot_motifs.jl](15_plot_motifs.jl) | ✅ Yes | ✅ Result files are available |
|  |  |  |
|  |  |  |
| **Preparation scripts** |  |  |
| (Usually no need to re-run as all outputs are on GitHub) |  |  |
| Any script in scripts/prep | ✅ Yes | ✅ Yes |
|  |  |  |
