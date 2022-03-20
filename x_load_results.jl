## Load data

# Define names for the 4 sampling options
options = ["mean", "mean_thr", "rand", "rand_thr"]
titles = ["Mean" "Mean > cutoff" "Rnd" "Rnd > cutoff"] # for plots later on

# Richness layers
S_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath("data", "results", "richness_$(opt).tif")
    S_all[opt] = geotiff(SimpleSDMPredictor, path)
end
S_all

# Species LCBD layers
lcbd_species_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath("data", "results", "lcbd_species_$(opt).tif")
    lcbd_species_all[opt] = geotiff(SimpleSDMPredictor, path)
end
lcbd_species_all

# Networks LCBD layers
lcbd_networks_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath("data", "results", "lcbd_networks_$(opt).tif")
    lcbd_networks_all[opt] = geotiff(SimpleSDMPredictor, path)
end
lcbd_networks_all

# Others
Sσ = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_uncertainty.tif"))
spatialrange = (left=-80.0, right=-50.0, bottom=45.0, top=65.)
reference_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; spatialrange...)