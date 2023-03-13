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

# Define a function to get the minimum non-zero value
function minimum_nonzero(x; dims=1)
    t = eltype(x) isa Union ? eltype(x).b : eltype(x)
    v = prevfloat(typemax(t))
    if ndims(x) > 1
        # For when we'll use the function on the ecoregion metaweb
        x_min = replace(minimum(replace(x, 0.0 => v); dims=dims), v => 0.0)
    else
        # For when we'll use the function on a summary ecoregion layer
        x_min = minimum(replace(x, 0.0 => v))
        x_min = isequal(x_min, v) ? 0.0 : x_min
    end
    return x_min
end

# Define the network measures to use
network_measures = ["Co", "L", "Lv", "Ld"]
network_fs = [connectance, links, links_var, linkage_density]
network_filename = ["connectance", "links_mean", "links_var", "links_density"]
summary_fs = [mean, median, maximum, minimum_nonzero]

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

#### Metaweb by ecoregion

#=

# Load layer with networks in each cell
include("09_get_network_measures.jl")

# Assemble ecoregion metaweb via the networks BitArray
function ecoregionalize(layer::T, networks::BitArray{4}, ecoregions_stack; fmeta=mean, fnet=links) where {T<: SimpleSDMResponse{UnipartiteProbabilisticNetwork{Float64, String}}}
    l_eco = similar(layer, Float32)
    @threads for e in ecoregions_stack
        e_inds = indexin(keys(e), keys(layer))
        networks_e = @view networks[e_inds, :, :, :];
        mat_e = dropdims(mean(networks_e; dims=4), dims=4)
        A_e = dropdims(fmeta(mat_e; dims=1), dims=1)
        network_e = UnipartiteProbabilisticNetwork(A_e, species(collect(layer)[1]))
        l_eco[keys(e)] = fill(fnet(network_e), length(keys(e)))
    end
    return l_eco
end

# Run for all options
ecometaweb_layers = Dict{String, SimpleSDMResponse}()
for o in opt
    ecometaweb_layers["$(o.m)_$(o.fm)"] = ecoregionalize(
        layer, networks, ecoregions_stack; fmeta=o.fm, fnet=o.fn
    )
end

# Export ecoregion metaweb results
for o in opt
    p = joinpath(ecoresults_path, "ecometaweb_$(o.m)_$(o.fm).tif")
    l = ecometaweb_layers["$(o.m)_$(o.fm)"]
    geotiff(p, l)
end

=#