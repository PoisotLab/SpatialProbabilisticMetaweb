#### Ecoregions

# CAN = true
include("A0_required.jl")

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    eco_path = joinpath("data", "input", "canada_ecoregions.tif");
    results_path = joinpath("data", "results");
else
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    eco_path = joinpath("data", "input", "quebec_ecoregions.tif")
    results_path = joinpath("xtras", "results");
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
    for e in ecoregions_stack
        l_eco[keys(e)] = fill(f(layer[keys(e)]), length(keys(e)))
    end
    if !keepzeros
        l_eco = replace(l_eco, 0.0 => nothing)
    end
    return l_eco
end

# Summarize by ecoregion
ecoregion_layers = Dict{String, SimpleSDMResponse}()
for m in network_measures
    ecoregion_layers[m] = ecoregionalize(local_layers[m], ecoregions_stack)
end

# Some variations
for f in [sum, median, maximum, minimum]
    ecoregion_layers["L_$(string(f))"] = ecoregionalize(
        local_layers["L"], ecoregions_stack; f=f
    )
end

#### Metaweb by ecoregion

# Load layer with networks in each cell
include("09_get_network_measures.jl")

# Assemble ecoregion metaweb via the networks BitArray
minimum_nonzero(x; dims=1) = replace(minimum(replace(x, 0.0 => 2.0); dims=dims), 2.0 => 0.0)
function ecoregionalize(layer::T, networks::BitArray{4}, ecoregions_stack; fmeta=mean, fnet=links) where {T<: SimpleSDMResponse{UnipartiteProbabilisticNetwork{Float64, String}}}
    # @assert fmeta in [mean, maximum, minimum, minimum_nonzero]
    l_eco = similar(layer, Float32)
    for e in ecoregions_stack
        e_inds = indexin(keys(e), keys(layer))
        networks_e = @view networks[e_inds, :, :, :];
        mat_e = dropdims(mean(networks_e; dims=4), dims=4)
        A_e = dropdims(fmeta(mat_e; dims=1), dims=1)
        replace!(A_e, 2.0 => 0.0)
        network_e = UnipartiteProbabilisticNetwork(A_e, species(layer[1]))
        l_eco[keys(e)] = fill(fnet(network_e), length(keys(e)))
    end
    return l_eco
end
ecometaweb_layers = Dict{String, SimpleSDMResponse}()
for (m,f) in zip(network_measures, network_fs)
    ecometaweb_layers[m] = ecoregionalize(
        layer, networks, ecoregions_stack; fmeta=mean, fnet=f
    )
end

# Test other options for metaweb assembly
for f in [mean, median, maximum, minimum]
    ecometaweb_layers["L_$(string(f))"] = ecoregionalize(
        layer, networks, ecoregions_stack; fmeta=f, fnet=links
    )
end