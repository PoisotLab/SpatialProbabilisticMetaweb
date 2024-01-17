#### Ecoregions

# CAN = true
include("A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    eco_path = joinpath("data", "input", "canada_ecoregions.tif");
    ref_path = joinpath("data", "input", "canada_ref_2.tif")
    results_path = joinpath("data", "results");
else
    eco_path = joinpath("data", "input", "quebec_ecoregions.tif")
    ref_path = joinpath("data", "input", "quebec_ref_10.tif")
    results_path = joinpath("xtras", "results");
end

# Load ecoregion objects and functions
include("scripts/lib/A5_ecoregions.jl")

## Basic summary statistics

# Define the network measures to use
measures = [
    "Co", "L", "Lv", "Ld", "S", "Sv", "LCBD_species", "LCBD_networks",
]
filenames = [
    "connectance", "links_mean", "links_var", "links_density",
    "richness_mean", "richness_uncertainty",
    "lcbd_species_mean", "lcbd_networks_mean",
]

# Load layers to summarize by ecoregion
local_layers = Dict{String, SimpleSDMPredictor}()
for (m, f) in zip(measures, filenames)
    local_layers[m] = read_geotiff(joinpath(results_path, "$f.tif"), SimpleSDMPredictor)
end

# Predefine set of options
opt = []
for m in measures, fs in summary_fs
    o = (m = m, fs = fs)
    push!(opt, o)
end
opt

# Summarize by ecoregion
ecoregion_layers = Dict{String, SimpleSDMResponse}()
@showprogress "Ecoregions:" for o in opt
    ecoregion_layers["$(o.m)_$(o.fs)"] = ecoregionalize(
        local_layers[o.m], ecoregions_stack; f=o.fs
    )
end
ecoregion_layers

# Rescale LCBD layers as relative values
for m in ["LCBD_species", "LCBD_networks"], f in summary_fs
    ecoregion_layers["$(m)_$f"] = ecoregion_layers["$(m)_$f"]/maximum(local_layers[m])
end

# Export layers
ecoresults_path = joinpath(results_path, "ecoregions");
isdir(ecoresults_path) || mkdir(ecoresults_path)
for (key, layer) in ecoregion_layers
    path = joinpath(ecoresults_path, "ecoregion_$key.tif")
    write_geotiff(path, layer)
end
