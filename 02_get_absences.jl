include("A0_required.jl")

# Option to use job array
# JOBARRAY = true

# Option to use WithinRadius for pseudoabsences
# WITHIN_RADIUS = true

# Option to run for CAN
# CAN = true
if (@isdefined CAN) && CAN == true
    occ_path = joinpath("data", "occurrences");
    pa_path = joinpath("data", "presence_absence");
    input_path = joinpath("data", "input");
    @info "Running for Canada at 2.5 arcmin resolution"
else
    occ_path = joinpath("xtras", "occurrences");
    pa_path = joinpath("xtras", "presence_absence");
    input_path = joinpath("xtras", "input");
    @info "Running for Quebec at 10 arcmin resolution"
end

reference_layer = read_geotiff(
    joinpath(input_path, "landcover_stack.tif"), SimpleSDMPredictor
)

occfiles = readdir(occ_path; join=true)
filter!(!contains("all_occurrences"), occfiles) # remove backup files for QC data
filter!(!contains(".gitkeep"), occfiles)

# Divide files if using job array
if (@isdefined JOBARRAY) && JOBARRAY == true
    _jobid = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
    _jobcount = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_COUNT", "1"))

    _ntot = length(occfiles)
    _nfiles = ceil(Int64, _ntot/_jobcount)

    # Balance larger files between jobs
    _gp = repeat(vcat(collect(1:_jobcount), collect(_jobcount:-1:1)), ceil(Int64, _nfiles/2))[1:_ntot]
    _gp_f = [sort(occfiles, by=filesize, rev=true)[findall(==(i), _gp)] for i in 1:_jobcount]

    jobfiles = _gp_f[_jobid]
else
    jobfiles = occfiles
end

ispath(pa_path) || mkpath(pa_path)

p = Progress(length(occfiles))

Threads.@threads for i in axes(jobfiles, 1)
    try
        pres = similar(reference_layer, Bool)
        df = DataFrame(CSV.File(jobfiles[i]; stringtype=String))
        for r in eachrow(df)
            if !isnothing(pres[r.longitude, r.latitude])
                pres[r.longitude, r.latitude] = true
            end
        end
        if (@isdefined WITHIN_RADIUS) && WITHIN_RADIUS == true
            absmask = pseudoabsencemask(WithinRadius, pres)
        else
            absmask = pseudoabsencemask(DistanceToEvent, pres) # * cellsize(pres)
        end
        Random.seed!(i)
        abs = backgroundpoints(absmask, sum(pres); replace=false)
        replace!(abs, false => nothing)
        replace!(pres, false => nothing)
        outfile = replace(jobfiles[i], "occurrences" => "presence_absence", ".csv" => ".tif")
        write_geotiff(outfile, [convert(Float32, pres), convert(Float32, abs)])
        if !(@isdefined quiet) || quiet == false
            # Print progress bar
            next!(p)
        end
    catch
        @warn "Error with $(jobfiles[i])"
    end
end
