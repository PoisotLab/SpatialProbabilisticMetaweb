#### Richness & LCBD analysis ####

# CAN = true
include("04_aggregate_sdms.jl");

# Load the corresponding results if dealing with CAN data or minimal example
if (@isdefined CAN) && CAN == true
    fit_path = joinpath("data", "input", "sdm_fit_results.csv")
    results_path = joinpath("data", "results")
else
    fit_path = joinpath("xtras", "input", "sdm_fit_results.csv")
    results_path = joinpath("xtras", "results")
end

## Richness

# We will create maps of species richness and standard deviation
# We do so by exploring 4 options to obtain the richness values

# Before everything let's extract the species list
spp = collect(keys(D))

# Option 1: Species richness from distribution means
Smeans = Dict{String, SimpleSDMPredictor}()
Svars = Dict{String, SimpleSDMPredictor}()
@threads for sp in spp
    Smeans[sp] = mean.(D[sp])
    Svars[sp] = var.(D[sp])
end
Sμ = reduce(+, values(Smeans))
Sv = reduce(+, values(Svars))

# Option 2: Species richness for random samples
Random.seed!(42)
Srands = Dict{String, SimpleSDMPredictor}()
@threads for sp in spp
    Srands[sp] = rand.(D[sp])
end
Sr = reduce(+, values(Srands))

# Option 3-4: Convert to presence absence data based on cutoff value
# Prepare cutoff values for all species
sdm_results = CSV.read(fit_path, DataFrame)
sdm_results.species = replace.(sdm_results.species, "_" => " ")
cutoffs = Dict{String, Float64}()
for r in eachrow(sdm_results)
    cutoffs[r.species] = r.cutoff
end

# Some species have no cutoffs, so let's add an impossible one to make everything work
missing_sp = setdiff(spp, sdm_results.species)
for m in missing_sp
    cutoffs[m] = 1.0
end
cutoffs

# We can now get the species richness for the thresholded distributions
Smeans_cut = Dict{String, SimpleSDMResponse}()
Srands_cut = Dict{String, SimpleSDMResponse}()
for sp in spp
    Smeans_cut[sp] = Smeans[sp] .> cutoffs[sp]
    Srands_cut[sp] = Srands[sp] .> cutoffs[sp]
end
Sμ_cut = reduce(+, convert.(Float64, values(Smeans_cut)))
Sr_cut = reduce(+, convert.(Float64, values(Srands_cut)))

# Export results
isdir(results_path) || mkpath(results_path)
write_geotiff(joinpath(results_path, "richness_mean.tif"), Sμ)
write_geotiff(joinpath(results_path, "richness_uncertainty.tif"), Sv)
write_geotiff(joinpath(results_path, "richness_rand.tif"), Sr)
write_geotiff(joinpath(results_path, "richness_mean_thr.tif"), Sμ_cut)
write_geotiff(joinpath(results_path, "richness_rand_thr.tif"), Sr_cut)

## LCBD values

# Get LCBD values for all 4 assembly options
isdir(results_path) || mkpath(results_path)
options = ["mean", "mean_thr", "rand", "rand_thr"]
lcbd_layers = Dict{String, SimpleSDMResponse}()
for (i, opt) in enumerate(options)
    S = values([Smeans, Srands, Smeans_cut, Srands_cut][i])
    @assert length(unique(length.(S))) == 1

    # Y matrix
    Y = reduce(hcat, values.(S))

    # Temporary fix for bug with negative values
    inds_neg = findall(<(0.0), Y) # 2 values only
    if length(inds_neg) > 0
        @info "$(length(inds_neg)) negative values were replaced by zero"
        @info Y[inds_neg]
        Y[inds_neg] .= 0.0
    end

    # LCBD
    lcbd_layers[opt] = similar(reference_layer)
    lcbd_layers[opt][keys(reference_layer)] = LCBD(hellinger(Y))[1]

    # Export results
    write_geotiff(joinpath(results_path, "lcbd_species_$opt.tif"), lcbd_layers[opt])
end
