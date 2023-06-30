#### Ecoregions

CAN = true
include("A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    eco_path = joinpath("data", "input", "canada_ecoregions_can.tif");
    ref_path = joinpath("data", "input", "canada_ref_2.tif")
    results_path = joinpath("data", "results");
else
    eco_path = joinpath("data", "input", "quebec_ecoregions_can.tif")
    ref_path = joinpath("data", "input", "quebec_ref_10.tif")
    results_path = joinpath("xtras", "results");
end

# Load ecoregions
ecoregions = read_geotiff(eco_path, SimpleSDMPredictor)
heatmap(ecoregions)

# Separate ecoregions in different layers
ecoregions_ids = unique(values(ecoregions))
ecoregions_stack = [convert(Float32, ecoregions .== id) for id in ecoregions_ids]
for e in ecoregions_stack
    replace!(e, 0.0 => nothing)
end

## Basic summary statistics

# Define the network measures to use
measures = ["Co", "L", "Lv", "Ld", "S", "Sv", "LCBD_species", "LCBD_networks"]
filenames = [
    "connectance", "links_mean", "links_var", "links_density",
    "richness_mean", "richness_uncertainty",
    "lcbd_species_mean", "lcbd_networks_mean"
]

# Define the summary functions we will use
quantile055(x) = quantile(x, 0.055)
quantile945(x) = quantile(x, 0.945)
iqr89(x) = quantile945(x) - quantile055(x)
summary_fs = [median, quantile055, quantile945, iqr89]

# Load layers to summarize by ecoregion
local_layers = Dict{String, SimpleSDMPredictor}()
for (m, f) in zip(measures, filenames)
    local_layers[m] = read_geotiff(joinpath(results_path, "$f.tif"), SimpleSDMPredictor)
end

# Define function
function ecoregionalize(layer, ecoregions_stack; f=median, keepzeros=true)
    l_eco = similar(layer)
    @threads for e in ecoregions_stack
        l_eco[keys(e)] = fill(f(layer[keys(e)]), length(e))
    end
    if !keepzeros
        replace!(l_eco, 0.0 => nothing)
    end
    return l_eco
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
for o in opt
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
for o in opt
    path = joinpath(ecoresults_path, "ecoregion_$(o.m)_$(o.fs).tif")
    layer = ecoregion_layers["$(o.m)_$(o.fs)"]
    write_geotiff(path, layer)
end
