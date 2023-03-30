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
    S_all[opt] = geotiff(SimpleSDMPredictor, path)
end
S_all
Sv = geotiff(SimpleSDMPredictor, joinpath(results_path, "richness_uncertainty.tif"))

# Species LCBD layers
lcbd_species_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath(results_path, "lcbd_species_$(opt).tif")
    lcbd_species_all[opt] = geotiff(SimpleSDMPredictor, path)
end
lcbd_species_all
# Temporarily replace values in the layer with NaNs
lcbd_species_all["rand_thr"] = replace(lcbd_species_all["rand_thr"], NaN => 0.0)

# Networks LCBD layers
lcbd_networks_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath(results_path, "lcbd_networks_$(opt).tif")
    if isfile(path)
        lcbd_networks_all[opt] = geotiff(SimpleSDMPredictor, path)
    end
end
lcbd_networks_all

# Set LCBD values to their relative maximum
for lcbd_set in [lcbd_species_all, lcbd_networks_all], opt in keys(lcbd_set)
    lcbd_set[opt] = lcbd_set[opt]/maximum(lcbd_set[opt])
end