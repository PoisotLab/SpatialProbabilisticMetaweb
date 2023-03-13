#### Ecoregions

# CAN = true
include("A0_required.jl")

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    eco_path = joinpath("data", "input", "canada_ecoregions.tif");
    results_path = joinpath("data", "results");
    ecoresults_path = joinpath("data", "ecoregions");
else
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    eco_path = joinpath("data", "input", "quebec_ecoregions.tif")
    results_path = joinpath("xtras", "results");
    ecoresults_path = joinpath("xtras", "ecoregions");
end

# Define reference layer
reference_layer = geotiff(SimpleSDMPredictor, ref_path)
spatialrange = boundingbox(reference_layer)

# Load ecoregions
ecoregions = geotiff(SimpleSDMPredictor, eco_path)
plot(ecoregions)

# Separate ecoregions in different layers
ecoregions_ids = unique(collect(ecoregions))
ecoregions_stack = [convert(Float32, broadcast(==(e), ecoregions)) for e in ecoregions_ids]
ecoregions_stack = [replace(e, 0.0 => nothing) for e in ecoregions_stack]

## Basic summary statistics

# Define the network measures to use
network_measures = ["Co", "L", "Lv", "Ld"]
network_fs = [connectance, links, links_var, linkage_density]
network_filename = ["connectance", "links_mean", "links_var", "links_density"]

# Define the summary functions we will use
quantile055(x) = quantile(x, 0.055)
quantile945(x) = quantile(x, 0.945)
iqr89(x) = quantile945(x) - quantile055(x)
summary_fs = [median, quantile055, quantile945, iqr89]

# Predefine set of options for threading
opt = []
for (m, fn) in zip(network_measures, network_fs), fm in summary_fs
    o = (m = m, fn = fn, fm = fm)
    push!(opt, o)
end

# Load layers to summarize by ecoregion
local_layers = Dict{String, SimpleSDMPredictor}()
for (m, f) in zip(network_measures, network_filename)
    local_layers[m] = geotiff(SimpleSDMPredictor, joinpath(results_path, "$f.tif"))
end
local_layers["S"] = geotiff(SimpleSDMPredictor, joinpath(results_path, "richness_mean.tif"))
local_layers["Sσ"] = geotiff(SimpleSDMPredictor, joinpath(results_path, "richness_uncertainty.tif"))

# Define function
function ecoregionalize(layer, ecoregions_stack; f=mean, keepzeros=true)
    l_eco = similar(layer)
    @threads for e in ecoregions_stack
        l_eco[keys(e)] = fill(f(layer[keys(e)]), length(keys(e)))
    end
    if !keepzeros
        l_eco = replace(l_eco, 0.0 => nothing)
    end
    return l_eco
end

# Summarize by ecoregion
ecoregion_layers = Dict{String, SimpleSDMResponse}()
for o in opt
    ecoregion_layers["$(o.m)_$(o.fm)"] = ecoregionalize(
        local_layers[o.m], ecoregions_stack; f=o.fm
    )
end

# Export layers
isdir(ecoresults_path) || mkdir(ecoresults_path)
for o in opt
    p = joinpath(ecoresults_path, "ecoregion_$(o.m)_$(o.fm).tif")
    l = ecoregion_layers["$(o.m)_$(o.fm)"]
    geotiff(p, l)
end
