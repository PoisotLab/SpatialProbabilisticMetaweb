using SimpleSDMLayers
using StatsPlots
using JSON
using EcologicalNetworks
using DataFrames
import CSV
using ProgressMeter
using Statistics
using Distributions
using Clustering
using MultivariateStats

include("A1_LCBD.jl")

spatialrange = (left=-80., right=-50., bottom=45., top=65.)
reference_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; spatialrange...)
plot(reference_layer)

map_files = joinpath.("sdms", readdir("sdms"))

μ = Dict{String,SimpleSDMPredictor}()
σ = Dict{String,SimpleSDMPredictor}()
for map_file in map_files
    if contains(map_file, "_error.tif")
        sp_name = sp_name = replace(replace(replace(map_file, "sdms/" => ""), "_" => " "), "error.tif" => "")[1:end-1]
        σ[sp_name] = clip(geotiff(SimpleSDMPredictor, map_file); spatialrange...)
    else
        sp_name = sp_name = replace(replace(replace(map_file, "sdms/" => ""), "_" => " "), "model.tif" => "")[1:end-1]
        μ[sp_name] = clip(geotiff(SimpleSDMPredictor, map_file); spatialrange...)
    end
end

# We need a few zero types for distributions, which will allow to use them in cell of layers
Base.zero(::Type{Normal{T}}) where T = Normal(zero(T), zero(T))
Base.zero(::Type{Bernoulli{T}}) where {T} = Bernoulli(zero(T))
Base.zero(::Type{Truncated{Normal{T}, Continuous, T}}) where {T} = Truncated(zero(Normal{T}), zero(T), one(T))

D = Dict{String, SimpleSDMResponse}()
for sp in keys(μ)
    _t = similar(μ[sp], Truncated{Normal{Float64}, Continuous, Float64})
    for site in keys(μ[sp])
        _t[site] = Truncated(Normal(Float64(μ[sp][site]), Float64(σ[sp][site])), 0.0, 1.0)
    end
    D[sp] = _t
    GC.gc()
end

# Map of species richness and standard deviation
Smeans = map(l -> broadcast(mean, l), collect(values(D)))
Sμ = reduce(+, Smeans)
Sσ = sqrt(reduce(+, map(l -> broadcast(var, l), collect(values(D)))))

plot(Sμ, title="Expected richness", c=:batlow)
plot(Sσ, title="Std. dev. of richness", c=:acton)

bivariate(Sμ, Sσ; quantiles=true, classes=3)

Y = zeros(Float64, (length(reference_layer), length(Smeans)))
for i in eachindex(Smeans)
    Y[:,i] = Smeans[i][keys(reference_layer)]
end

lcbd_species = similar(reference_layer)
lcbd_species[keys(reference_layer)] = LCBD(hellinger(Y))[1]
plot(lcbd_species)

# Download and parse the metaweb
mw_output = DataFrame(CSV.File(download("https://raw.githubusercontent.com/PoisotLab/MetawebTransferLearning/main/artifacts/canadian_thresholded.csv"); stringtype=String))
sp = collect(keys(μ))
S = fill(0.0, length(sp), length(sp))
P = UnipartiteProbabilisticNetwork(S, sp)
for r in eachrow(mw_output)
    from = findfirst(isequal(r.from), sp)
    to = findfirst(isequal(r.to), sp)
    P[from, to] = r.score
end

A = zeros(Float64, richness(P), richness(P))
for (i,si) in enumerate(species(P))
    for (j,sj) in enumerate(species(P))
        A[i,j] = P[si,sj]
    end
end

# Get a preliminary map
sites = keys(reference_layer)
networks = zeros(Bool, length(sites), size(P)..., 10)
p = Progress(length(sites))
Threads.@threads for i in eachindex(sites)
    site = sites[i]
    s = [D[s][site] for s in species(P)]
    pcooc = mean.(s) .* mean.(s)'
    # 1 [sdm - | classifier -] --> avg * avg'
    # 2 [sdm - | classifier +] --> (avg > thr) * (avg > thr)'
    # 3 [sdm + | classifier -] --> rnd * rnd'
    # 4 [sdm + | classifier +] --> (rnd > thr) * (rnd > thr)'
    for j in 1:size(networks, 4)
        networks[i, :, :, j] .= adjacency(rand(UnipartiteProbabilisticNetwork(pcooc .* A, species(P))))
    end
    next!(p)
end

valued_interactions = findall(!iszero, sum(networks; dims=(1,4))[1,:,:])
by_site = sum(networks; dims=(4))

Z = zeros(Int64, (length(reference_layer), length(valued_interactions)))
for i in 1:length(sites)
    Z[i,:] = by_site[i,valued_interactions]
end

lcbd_networks = similar(reference_layer)
lcbd_networks[keys(reference_layer)] = LCBD(hellinger(Z))[1]

plot(
    plot(lcbd_species, leg=false),
    plot(lcbd_networks, leg=false)
)

p0 = colorant"#e8e8e8"
bv_pal_1 = (p0=p0, p1=colorant"#64acbe", p2=colorant"#c85a5a")
bv_pal_2 = (p0=p0, p1=colorant"#73ae80", p2=colorant"#6c83b5")
bv_pal_3 = (p0=p0, p1=colorant"#9972af", p2=colorant"#c8b35a")
bv_pal_4 = (p0=p0, p1=colorant"#be64ac", p2=colorant"#5ac8c8")

bivariate(lcbd_networks, lcbd_species; quantiles=true, bv_pal_1..., classes=5)
plot(rescale(lcbd_networks, collect(0.0:0.05:1.0)); c=cgrad([p0, colorant"#be64ac"]))
plot(rescale(lcbd_species, collect(0.0:0.05:1.0)); c=cgrad([p0, colorant"#5ac8c8"]))

histogram2d(
    rescale(lcbd_networks, collect(0.0:0.05:1.0)),
    rescale(lcbd_species, collect(0.0:0.05:1.0));
    bins=20
)
xaxis!((0,1))
yaxis!((0,1))