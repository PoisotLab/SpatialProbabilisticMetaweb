#### Main analysis script ####

# Load required packages and define custom functions
include("A0_required.jl");

# Options to control outputs
# JOBARRAY = true # whether to run as a job array or not
# CAN = true # to run script for Canada, which is resource-intensive
# quiet = true # whether to disable progress bar (mostly for clusters)

## Preparation scripts

# These scripts are not part of the analysis pipeline but produce required
# elements which are version-controlled and available on the repo (e.g. in the
# data/input/ folder).

include("scripts/prep/P1_get_canada_shapefile.jl");
include("scripts/prep/P2_get_climate.jl");
include("scripts/prep/P3_get_landcover.jl");
include("scripts/prep/P4_get_ecoregions.jl");
include("scripts/prep/P5_reconcile_metaweb.jl");

## Part I - Prepare data ####

# Get occurrences from GBIF
include("01_get_occurrences.jl");

# Create pseudo-absence
# WITHIN_RADIUS = true # to use WithinRadius option (needs to run on Narval for Canada)
include("02_get_absences.jl");

# Train SDMs and predict species distributions
include("03_generate_sdms.jl");

## Part II - Assemble results ####

# These scripts are resource-intensive at the Canada scale and need to run on
# the clusters.

# These scripts are run by later scripts to assemble the species distribution
# and network results (which are too large to be exported). There's usually no
# need to run them alone.
include("04_aggregate_sdms.jl"); # faster, not memory-intensive
include("05_assemble_networks.jl"); # longer, memory-intensive, requires script 04

# Get the species LCBD
include("06_get_species_lcbd.jl"); # requires script 04

# Get the network LCBD
include("07_get_network_lcbd.jl"); # requires scripts 04 and 05

# Get the network measures
include("09_get_network_measures.jl"); # requires scripts 04 and 05

# Get the ecoregions measures (can run locally)
include("11_get_ecoregions_measures.jl"); # requires the results of 06, 07, 09

## Part III - Plot results ####

# These scripts use exported layers and can run locally

# Set SAVE to true to export figures with CairoMakie
SAVE = true

# Plot LCBD results (species and networks)
include("08_plot_lcbd_results.jl");

# Plot network measures
include("10_plot_network_measures.jl");

# Plot ecoregion results
include("12_plot_ecoregions.jl");

## Additional analyses

# These has-beens were actually really cool in their own time, you know. We
# can't throw away our past just like that.

# include("xtras/archive/x_get_network_extras.jl");
# include("xtras/archive/x_get_ecoregions_metaweb.jl");
# include("xtras/archive/x_working_group.jl");
