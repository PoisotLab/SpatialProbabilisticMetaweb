#### Richness & LCBD analysis ####

# QC = true
include("04_aggregate_sdms.jl")

# Load the previous sdm results if dealing with QC data
if (@isdefined QC) && QC == true
    fit_path = joinpath("xtras", "input", "sdm_fit_results.csv")
    results_path = joinpath("xtras", "results")
else
    fit_path = joinpath("data", "input", "sdm_fit_results.csv")
    results_path = joinpath("data", "results")
end

## Richness

# We will create maps of species richness and standard deviation
# We do so by exploring 4 options to obtain the richness values

# Option 1: Species richness from distribution means
Smeans = map(l -> broadcast(mean, l), collect(values(D)))
Sμ = reduce(+, Smeans)
Svars = map(l -> broadcast(var, l), collect(values(D)))
Sσ = sqrt(reduce(+, Svars))

# Option 2: Species richness for random samples
Random.seed!(42)
Srands = map(l -> broadcast(rand, l), collect(values(D)))
Sr = reduce(+, Srands)

# Option 3-4: Convert to presence absence data based on cutoff value
# Prepare cutoff values for all species
sdm_results = CSV.read(fit_path, DataFrame)
sdm_results.species = replace.(sdm_results.species, "_" => " ")
cutoffs = Dict{String, Float64}()
for r in eachrow(sdm_results)
    cutoffs[r.species] = r.cutoff
end

# Some species have no cutoffs, so let's add an impossible one to make everything work
missing_sp = setdiff(keys(D), sdm_results.species)
for m in missing_sp
    cutoffs[m] = 1.0
end
cutoffs

# We can now get the species richness for the thresholded distributions
spp = collect(keys(D))
Smeans_cut = [broadcast(>(cutoffs[sp]), S) for (sp, S) in zip(spp, Smeans)]
Srands_cut = [broadcast(>(cutoffs[sp]), S) for (sp, S) in zip(spp, Srands)]
Sμ_cut, Sr_cut = [reduce(+, convert.(Float32, S)) for S in (Smeans_cut, Srands_cut)]

# Export results
isdir(results_path) || mkpath(results_path)
geotiff(joinpath(results_path, "richness_mean.tif"), Sμ)
geotiff(joinpath(results_path, "richness_uncertainty.tif"), Sσ)
geotiff(joinpath(results_path, "richness_rand.tif"), Sr)
geotiff(joinpath(results_path, "richness_mean_thr.tif"), Sμ_cut)
geotiff(joinpath(results_path, "richness_rand_thr.tif"), Sr_cut)

## LCBD values

# Get LCBD values for all 4 assembly options
lcbd_layers = []
for S in (Smeans, Srands, Smeans_cut, Srands_cut)
    # Y matrix
    Y = zeros(Float64, (length(reference_layer), length(S)))
    for i in eachindex(S)
        Y[:,i] = S[i][keys(reference_layer)]
    end

    # LCBD
    lcbd_species = similar(reference_layer)
    lcbd_species[keys(reference_layer)] = LCBD(hellinger(Y))[1]
    push!(lcbd_layers, lcbd_species)
end

# Export results
geotiff(joinpath(results_path, "lcbd_species_mean.tif"), lcbd_layers[1])
geotiff(joinpath(results_path, "lcbd_species_rand.tif"), lcbd_layers[2])
geotiff(joinpath(results_path, "lcbd_species_mean_thr.tif"), lcbd_layers[3])
geotiff(joinpath(results_path, "lcbd_species_rand_thr.tif"), lcbd_layers[4])