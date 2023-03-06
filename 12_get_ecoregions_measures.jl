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

# Load layers to summarize by ecoregion
Co = geotiff(SimpleSDMPredictor, joinpath(results_path, "connectance.tif"))
L = geotiff(SimpleSDMPredictor, joinpath(results_path, "links_mean.tif"))
Lv = geotiff(SimpleSDMPredictor, joinpath(results_path, "links_var.tif"))
Ld = geotiff(SimpleSDMPredictor, joinpath(results_path, "links_density.tif"))
S = geotiff(SimpleSDMPredictor, joinpath(results_path, "richness_mean.tif"))
Sσ = geotiff(SimpleSDMPredictor, joinpath(results_path, "richness_uncertainty.tif"))

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
Co_eco = ecoregionalize(Co, ecoregions_stack)
L_eco = ecoregionalize(L, ecoregions_stack)
Lv_eco = ecoregionalize(Lv, ecoregions_stack)
Ld_eco = ecoregionalize(Ld, ecoregions_stack)
S_eco = ecoregionalize(S, ecoregions_stack)
Sσ_eco = ecoregionalize(Sσ, ecoregions_stack)

#### Metaweb by ecoregion

# Load layer with networks in each cell
include("09_get_network_measures.jl")

# Assemble ecoregion metaweb via the networks BitArray
function ecoregionalize(layer::T, networks::BitArray{4}, ecoregions_stack; fsite=mean, fnet=links) where {T<: SimpleSDMResponse{UnipartiteProbabilisticNetwork{Float64, String}}}
    @assert fsite in [mean, max, min]
    l_eco = similar(layer, Float32)
    for e in ecoregions_stack
        e_inds = indexin(keys(e), keys(layer))
        networks_e = @view networks[e_inds, :, :, :];
        mat_e = dropdims(mean(networks_e; dims=4), dims=4)
        if fsite == mean
            A_e = dropdims(fsite(mat_e; dims=1), dims=1)
        elseif fsite in [max, min]
            A_e = dropdims(reduce(fsite, mat_e; dims=1), dims=1)
        end
        network_e = UnipartiteProbabilisticNetwork(A_e, species(layer[1]))
        l_eco[keys(e)] = fill(fnet(network_e), length(keys(e)))
    end
    return l_eco
end
Co_meta_eco = ecoregionalize(layer, networks, ecoregions_stack; fsite=mean, fnet=connectance)
L_meta_eco = ecoregionalize(layer, networks, ecoregions_stack; fsite=mean, fnet=links)
Lv_meta_eco = ecoregionalize(layer, networks, ecoregions_stack; fsite=mean, fnet=links_var)
Ld_meta_eco = ecoregionalize(layer, networks, ecoregions_stack; fsite=mean, fnet=linkage_density)

# Test other options for metaweb assembly
L_meta_eco_mean = ecoregionalize(layer, networks, ecoregions_stack; fsite=mean, fnet=links)
L_meta_eco_max = ecoregionalize(layer, networks, ecoregions_stack; fsite=max, fnet=links)
L_meta_eco_min = ecoregionalize(layer, networks, ecoregions_stack; fsite=min, fnet=links)
