#### Main analysis script ####

# Load required packages and define custom functions
include("A0_required.jl");

# Options to control outputs
# JOBARRAY = true # whether to run as a job array or not
# CAN = true # to run script for Canada, which is resource-intensive
# quiet = true # whether to disable progress bar (mostly for clusters)

## Part I - Prepare data ####

# Get occurrences from GBIF
@info "Running 01_get_occurrences.jl"
include("01_get_occurrences.jl");

# Create pseudo-absence
# WITHIN_RADIUS = true # to use WithinRadius option (needs to run on Narval for Canada)
@info "Running 02_get_absences.jl"
include("02_get_absences.jl");

# Train SDMs and predict species distributions
@info "Running 03_generate_sdms.jl"
include("03_generate_sdms.jl");

## Part II - Assemble results ####

# These scripts are resource-intensive at the Canada scale and need to run on
# the clusters.

# These scripts are run by following scripts to assemble the species
# distribution and network results (which are too large to be exported). There's
# usually no need to run them alone.
@info "Running scripts 04 & 05"
include("04_aggregate_sdms.jl"); # faster, not memory-intensive
include("05_assemble_networks.jl"); # longer, memory-intensive, requires script 04

# Get the species LCBD
@info "Running 06_get_species_lcbd.jl"
include("06_get_species_lcbd.jl"); # requires script 04

# Get the network LCBD
@info "Running script 07_get_network_lcbd.jl"
include("07_get_network_lcbd.jl"); # requires scripts 04 and 05

# Get the network measures
@info "Running 09_get_network_measures.jl"
include("09_get_network_measures.jl"); # requires scripts 04 and 05

# Get the ecoregions measures (can run locally)
@info "Running 11_get_ecoregions_measures.jl"
include("11_get_ecoregions_measures.jl"); # requires the results of 06, 07, 09

## Part III - Plot results ####

# These scripts use exported layers and can run locally

# Set SAVE to true to export figures with CairoMakie
# SAVE = true

# Plot LCBD results (species and networks)
@info "Running 08_plot_lcbd_results.jl"
@time include("08_plot_lcbd_results.jl");

# Plot network measures
@info "Running 10_plot_network_measures.jl"
include("10_plot_network_measures.jl");

# Plot ecoregion results
@info "Running 12_plot_ecoregions.jl"
include("12_plot_ecoregions.jl");

## Part IV - Motifs analysis ####

# This analysis is especially resource-intensive so we consider it separately

# Get motifs
MOTIF = :S4 # Need to choose a motif (S1, S2, S4, S5) and re-run separately for each
@info "Running 13_get_motifs.jl"
include("13_get_motifs.jl")

# Assemble motif result layers
@info "Running 14_assemble_motifs.jl"
include("14_assemble_motifs.jl")

# Plot motif results
@info "Running 15_plot_motifs.jl"
include("15_plot_motifs.jl")

## Preparation scripts ####

# These scripts are not part of the analysis pipeline but produce required
# elements which are version-controlled and available on the repo (e.g. in the
# data/input/ folder).

#=
include("scripts/prep/P1_get_canada_shapefile.jl");
include("scripts/prep/P2_get_climate.jl");
include("scripts/prep/P3_get_landcover.jl");
include("scripts/prep/P4_get_ecoregions.jl");
include("scripts/prep/P5_reconcile_metaweb.jl");
=#