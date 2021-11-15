using SimpleSDMLayers
using StatsPlots
using EcologicalNetworks
using DataFrames
import CSV
using ProgressMeter
using Statistics
using Distributions
using Random

include("A1_LCBD.jl")

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

# Map of species richness and standard deviation
Smeans = map(l -> broadcast(mean, l), collect(values(D)))
Sμ = reduce(+, Smeans)
Svars = map(l -> broadcast(var, l), collect(values(D)))
Sσ = sqrt(reduce(+, Svars))

# Univariate maps
p1 = plot(Sμ, title="Expected richness", c=:batlow)
p2 = plot(Sσ, title="Std. dev. of richness", c=:acton)

# Bivariate map
p3 = bivariate(Sμ, Sσ; quantiles=true, classes=3, xlab="Longitude", ylab="Latitude")
bivariatelegend!(
    Sμ,
    Sσ;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Expected richness",
    ylab="Std. dev. of richness",
    guidefontsize=7,
)

plot(p1, p2, p3, layout=(3,1), size=(600,1200))

# Y matrix
Y = zeros(Float64, (length(reference_layer), length(Smeans)))
for i in eachindex(Smeans)
    Y[:,i] = Smeans[i][keys(reference_layer)]
end

# LCBD
lcbd_species = similar(reference_layer)
lcbd_species[keys(reference_layer)] = LCBD(hellinger(Y))[1]
plot(lcbd_species)
