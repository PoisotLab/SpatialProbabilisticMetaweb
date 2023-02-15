## Load data

# CAN = true

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

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

# Species LCBD layers
lcbd_species_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath(results_path, "lcbd_species_$(opt).tif")
    lcbd_species_all[opt] = geotiff(SimpleSDMPredictor, path)
end
lcbd_species_all

# Networks LCBD layers
lcbd_networks_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath(results_path, "lcbd_networks_$(opt).tif")
    lcbd_networks_all[opt] = geotiff(SimpleSDMPredictor, path)
end
lcbd_networks_all