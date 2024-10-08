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
reference_layer = read_geotiff(ref_path, SimpleSDMResponse)
spatialrange = boundingbox(reference_layer)

# Divide sites if using job array
if (@isdefined JOBARRAY) && JOBARRAY == true
    # Get infos on the job array
    _jobid = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
    if @isdefined(TOTAL_JOBS) && !isnothing(TOTAL_JOBS)
        _jobcount = TOTAL_JOBS
    else
        _jobcont = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_COUNT", "10"))
    end

    # Assign sites to all tasks
    _ntot = length(reference_layer)
    _nsites = ceil(Int64, _ntot/_jobcount)
    _jobsites = repeat(1:_jobcount; inner=_nsites)

    # Set values to nothing for sites that are NOT for the current task
    _outsites = findall(!=(_jobid), _jobsites)
    filter!(<=(_ntot), _outsites) # make sure we don't go over the total
    _outkeys = keys(reference_layer)[_outsites]
    reference_layer[_outkeys] = fill(nothing, length(_outkeys))
end

# Select files to load
sdm_files = readdir(sdm_path; join=true)
filter!(!contains(".gitkeep"), sdm_files)
if iszero(length(sdm_files))
    prev = "03_generate_sdms.jl"
    @warn "Missing necessary files. Attempting to re-run previous script $prev"
    include(prev)
    sdm_files = readdir(sdm_path; join=true)
    filter!(!contains(".gitkeep"), sdm_files)
end

# Load predictions mean & variance layers
μ = Dict{String,SimpleSDMPredictor}()
σ = Dict{String,SimpleSDMPredictor}()
p = Progress(length(sdm_files), "Loading SDMs")
for map_file in sdm_files
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
p = Progress(length(μ), "Assembling distributions")
for sp in String.(keys(μ))
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
