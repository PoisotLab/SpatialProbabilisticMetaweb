#### Richness & LCBD analysis ####

include("04_aggregate_sdms.jl")

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
sdm_results = CSV.read(joinpath("data", "input", "sdm_fit_results.csv"), DataFrame)
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
geotiff(joinpath("data", "results", "richness_mean.tif"), Sμ)
geotiff(joinpath("data", "results", "richness_uncertainty.tif"), Sσ)
geotiff(joinpath("data", "results", "richness_rand.tif"), Sr)
geotiff(joinpath("data", "results", "richness_mean_thr.tif"), Sμ_cut)
geotiff(joinpath("data", "results", "richness_rand_thr.tif"), Sr_cut)

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

# Select the first one as the main result
lcbd_species = lcbd_layers[1]

# Export results
geotiff(joinpath("data", "results", "lcbd_species_mean.tif"), lcbd_layers[1])
geotiff(joinpath("data", "results", "lcbd_species_rand.tif"), lcbd_layers[2])
geotiff(joinpath("data", "results", "lcbd_species_mean_thr.tif"), lcbd_layers[3])
geotiff(joinpath("data", "results", "lcbd_species_rand_thr.tif"), lcbd_layers[4])

## Plot results

# All richness options
clim1 = mapreduce(minimum, min, [Sμ, Sr, Sμ_cut, Sr_cut])
clim2 = mapreduce(maximum, max, [Sμ, Sr, Sμ_cut, Sr_cut])
lims = (clim1, clim2)
plot(
    plot(Sμ; c=:cividis, title="Sμ", clim=lims),
    plot(Sr; c=:cividis, title="Sr", clim=lims),
    plot(Sμ_cut; c=:cividis, title="Sμ_cut", clim=lims),
    plot(Sr_cut; c=:cividis, title="Sr_cut", clim=lims);
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "richness_all.png"))

# Univariate richness maps
plot(
    plot(Sμ, title="Expected richness", c=cgrad([p0, bv_pal_2[2]])),
    plot(Sσ, title="Std. dev. of richness", c=cgrad([p0, bv_pal_2[3]]));
    layout=(2,1),
    size=(600, 600)
)
savefig(joinpath("figures", "richness_two-panels.png"))

# Bivariate richness map
bivariate(Sμ, Sσ; quantiles=true, classes=3, xlab="Longitude", ylab="Latitude", bv_pal_2...)
bivariatelegend!(
    Sμ,
    Sσ;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Expected richness",
    ylab="Std. dev. of richness",
    guidefontsize=7,
    bv_pal_2...
)
plot!(title=["Richness & uncertainty bivariate" ""])
savefig(joinpath("figures", "richness_bivariate.png"))

# Single LCBD map
plot(lcbd_species; c=:viridis, title="Species LCBD")
savefig(joinpath("figures", "lcbd_species.png"))

# All LCBD options
plot(
    plot(lcbd_layers[1]; c=:viridis, title="LCBD means"),
    plot(lcbd_layers[2]; c=:viridis, title="LCBD rands"),
    plot(lcbd_layers[3]; c=:viridis, title="LCBD means cut"),
    plot(lcbd_layers[4]; c=:viridis, title="LCBD rands cut"),
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "lcbd_species_all.png"))