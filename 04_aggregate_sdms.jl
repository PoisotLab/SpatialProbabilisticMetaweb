#### Aggregate SDMs ####

include("A0_required.jl")

## Probabilistic distributions

# Define reference layer
spatialrange = (left=-80., right=-50., bottom=45., top=65.)
reference_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; spatialrange...)

# Select files to load
map_files = joinpath.("sdms", readdir("sdms"))

# Load predictions mean & variance layers
μ = Dict{String,SimpleSDMPredictor}()
σ = Dict{String,SimpleSDMPredictor}()
for map_file in map_files
    sp_name = @chain map_file begin
        replace("sdms/" => "")
        replace("_" => " ")
        replace(" error.tif" => "")
        replace(" model.tif" => "")
    end
    if contains(map_file, "error.tif")
        σ[sp_name] = geotiff(SimpleSDMPredictor, map_file; spatialrange...)
    else
        μ[sp_name] = geotiff(SimpleSDMPredictor, map_file; spatialrange...)
    end
end

# We need a few zero types for distributions, which will allow to use them in cell of layers
Base.zero(::Type{Normal{T}}) where T = Normal(zero(T), zero(T))
Base.zero(::Type{Bernoulli{T}}) where {T} = Bernoulli(zero(T))
Base.zero(::Type{Truncated{Normal{T}, Continuous, T}}) where {T} = Truncated(zero(Normal{T}), zero(T), one(T))

# Create layers of Truncated Normal distributions given the mean & variance
D = Dict{String, SimpleSDMResponse}()
for sp in String.(keys(μ))
    _t = similar(μ[sp], Truncated{Normal{Float64}, Continuous, Float64})
    for site in keys(μ[sp])
        _t[site] = Truncated(Normal(Float64(μ[sp][site]), Float64(σ[sp][site])), 0.0, 1.0)
    end
    D[sp] = _t
end
GC.gc()
