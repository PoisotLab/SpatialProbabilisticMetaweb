## Load LCBD results

# CAN = true

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end;

## Load data

# Define names for the 4 sampling options
options = ["mean", "mean_thr", "rand", "rand_thr"]
titles = ["Mean" "Mean > cutoff" "Rnd" "Rnd > cutoff"] # for plots later on

# Richness layers
S_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath(results_path, "richness_$(opt).tif")
    S_all[opt] = read_geotiff(path, SimpleSDMPredictor)
end
S_all
Sv = read_geotiff(joinpath(results_path, "richness_uncertainty.tif"), SimpleSDMPredictor)

# Species LCBD layers
lcbd_species_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath(results_path, "lcbd_species_$(opt).tif")
    lcbd_species_all[opt] = read_geotiff(path, SimpleSDMPredictor)
end
lcbd_species_all
# Temporarily replace values in the layer with NaNs
lcbd_species_all["rand_thr"] = replace(lcbd_species_all["rand_thr"], NaN => 0.0)

# Networks LCBD layers
lcbd_networks_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath(results_path, "lcbd_networks_$(opt).tif")
    if isfile(path)
        lcbd_networks_all[opt] = read_geotiff(path, SimpleSDMPredictor)
    end
end
lcbd_networks_all

# Set LCBD values to their relative maximum
for lcbd_set in [lcbd_species_all, lcbd_networks_all], opt in keys(lcbd_set)
    lcbd_set[opt] = lcbd_set[opt]/maximum(lcbd_set[opt])
end

# We need to fix an issue with the network LCBN layers before we compare with species LCBD
# Some sites had no links, so their LCBD values was set to nothing to avoid NaNs everywhere
# Now we'll also set them to NaN for species LCBD to compare the rest of the two layers
if length(lcbd_species_all["mean"]) > length(lcbd_networks_all["mean"])
    _ndiff = length(lcbd_species_all["mean"]) - length(lcbd_networks_all["mean"])
    @info "Creating a species LCBD layers without $_ndiff sites with missing network LCBD values"
    _nan_sites = setdiff(keys(lcbd_species_all["mean"]), keys(lcbd_networks_all["mean"]))
    lcbd_species_nan = convert(SimpleSDMResponse, lcbd_species_all["mean"])
    lcbd_species_nan[_nan_sites] = fill(nothing, length(_nan_sites))
else
    lcbd_species_nan = lcbd_species_all["mean"]
end;