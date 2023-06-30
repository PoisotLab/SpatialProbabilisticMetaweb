#### Aggregate SDMs ####

include("A0_required.jl");

# Option to run for CAN
# CAN = true
if (@isdefined CAN) && CAN == true
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    sdm_path = joinpath("data", "sdms");
    @info "Running for Canada at 2.5 arcmin resolution"
else
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    sdm_path = joinpath("xtras", "sdms");
    @info "Running for Quebec at 10 arcmin resolution"
end

## Probabilistic distributions

# Define reference layer
reference_layer = read_geotiff(ref_path, SimpleSDMPredictor)
spatialrange = boundingbox(reference_layer)

# Select files to load
map_files = readdir(sdm_path; join=true)

# Load predictions mean & variance layers
μ = Dict{String,SimpleSDMPredictor}()
σ = Dict{String,SimpleSDMPredictor}()
p = Progress(length(map_files))
@threads for map_file in map_files
    sp_name = @chain map_file begin
        basename
        replace("_error.tif" => "")
        replace("_model.tif" => "")
        replace("_" => " ")
    end
    if contains(map_file, "error.tif")
        σ[sp_name] = read_geotiff(map_file, SimpleSDMPredictor; spatialrange...)
        σ[sp_name] = mask(reference_layer, σ[sp_name])
    else
        μ[sp_name] = read_geotiff(map_file, SimpleSDMPredictor; spatialrange...)
        μ[sp_name] = mask(reference_layer, μ[sp_name])
    end
    if !(@isdefined quiet) || quiet == false
        # Print progress bar
        next!(p)
    end
end

# We need a few zero types for distributions, which will allow to use them in cell of layers
Base.zero(::Type{Normal{T}}) where T = Normal(zero(T), zero(T))
Base.zero(::Type{Bernoulli{T}}) where {T} = Bernoulli(zero(T))
Base.zero(::Type{Truncated{Normal{T}, Continuous, T, T, T}}) where {T} = Truncated(zero(Normal{T}), zero(T), one(T))

# Create layers of Truncated Normal distributions given the mean & variance
D = Dict{String, SimpleSDMResponse}()
p = Progress(length(μ))
@threads for sp in String.(keys(μ))
    _t = similar(μ[sp], Truncated{Normal{Float64}, Continuous, Float64, Float64, Float64})
    for site in keys(μ[sp])
        _t[site] = Truncated(Normal(Float64(μ[sp][site]), Float64(σ[sp][site])), 0.0, 1.0)
    end
    D[sp] = _t
    if !(@isdefined quiet) || quiet == false
        # Print progress bar
        next!(p)
    end
end
GC.gc()
