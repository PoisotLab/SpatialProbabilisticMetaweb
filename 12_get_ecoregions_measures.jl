#### Ecoregions

CAN = true
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

# Plot results
plot(Co_eco; title="Co")
plot(L_eco; title="L")
plot(Lv_eco; title="Lv")
plot(Ld_eco; title="Ld")
plot(S_eco; title="S")
plot(Sσ_eco; title="Sσ")

# Some variations
plot(ecoregionalize(S, ecoregions_stack; f=sum))
plot(ecoregionalize(S, ecoregions_stack; f=maximum))
plot(ecoregionalize(S, ecoregions_stack; f=minimum))