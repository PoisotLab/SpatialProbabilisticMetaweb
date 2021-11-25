using SimpleSDMLayers
using StatsPlots
using EcologicalNetworks
using DataFrames
import CSV
using ProgressMeter
using Statistics
using Distributions
using Random

default(; dpi=200)

include("A1_LCBD.jl")

## Probabilistic distributions

# Define reference layer
spatialrange = (left=-80., right=-50., bottom=45., top=65.)
reference_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; spatialrange...)
plot(reference_layer)

map_files = joinpath.("sdms", readdir("sdms"))

# Load predictions mean & variance layers
μ = Dict{String,SimpleSDMPredictor}()
σ = Dict{String,SimpleSDMPredictor}()
Threads.@threads for map_file in map_files
    if contains(map_file, "_error.tif")
        sp_name = replace(replace(replace(map_file, "sdms/" => ""), "_" => " "), "error.tif" => "")[1:end-1]
        σ[sp_name] = geotiff(SimpleSDMPredictor, map_file; spatialrange...)
    else
        sp_name = replace(replace(replace(map_file, "sdms/" => ""), "_" => " "), "model.tif" => "")[1:end-1]
        μ[sp_name] = geotiff(SimpleSDMPredictor, map_file; spatialrange...)
    end
end

# We need a few zero types for distributions, which will allow to use them in cell of layers
Base.zero(::Type{Normal{T}}) where T = Normal(zero(T), zero(T))
Base.zero(::Type{Bernoulli{T}}) where {T} = Bernoulli(zero(T))
Base.zero(::Type{Truncated{Normal{T}, Continuous, T}}) where {T} = Truncated(zero(Normal{T}), zero(T), one(T))

# Create layers of Truncated Normal distributions given the mean & variance
D = Dict{String, SimpleSDMResponse}()
Threads.@threads for sp in String.(keys(μ))
    _t = similar(μ[sp], Truncated{Normal{Float64}, Continuous, Float64})
    for site in keys(μ[sp])
        _t[site] = Truncated(Normal(Float64(μ[sp][site]), Float64(σ[sp][site])), 0.0, 1.0)
    end
    D[sp] = _t
end
GC.gc()

## Richness

# Map of species richness and standard deviation
Smeans = map(l -> broadcast(mean, l), collect(values(D)))
Sμ = reduce(+, Smeans)
Svars = map(l -> broadcast(var, l), collect(values(D)))
Sσ = sqrt(reduce(+, Svars))

# Species richness for random samples
Random.seed!(42)
Srands = map(l -> broadcast(rand, l), collect(values(D)))
Sr = reduce(+, Srands)

# Species richness for thresholded distributions
spp = collect(keys(D))
Smeans_cut = [broadcast(>(cutoffs[sp]), S) for (sp, S) in zip(spp, Smeans)]
Srands_cut = [broadcast(>(cutoffs[sp]), S) for (sp, S) in zip(spp, Srands)]
Sμ_cut, Sr_cut = [reduce(+, convert.(Float32, S)) for S in (Smeans_cut, Srands_cut)]

# Plot all options
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

# Prepare colors
p0 = colorant"#e8e8e8"
bv_pal_2 = (p0=p0, p1=colorant"#73ae80", p2=colorant"#6c83b5")

# Univariate maps
plot(
    plot(Sμ, title="Expected richness", c=cgrad([p0, bv_pal_2[2]])),
    plot(Sσ, title="Std. dev. of richness", c=cgrad([p0, bv_pal_2[3]]));
    layout=(2,1),
    size=(600, 600)
)
savefig("figures/richness_two-panels.png")

# Bivariate map
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

## LCBD values

# Y matrix
Y = zeros(Float64, (length(reference_layer), length(Smeans)))
for i in eachindex(Smeans)
    Y[:,i] = Smeans[i][keys(reference_layer)]
end

# LCBD
lcbd_species = similar(reference_layer)
lcbd_species[keys(reference_layer)] = LCBD(hellinger(Y))[1]

# Plot LCBD
plot(lcbd_species; c=:viridis, title="Species LCBD")
savefig(joinpath("figures", "lcbd_species.png"))