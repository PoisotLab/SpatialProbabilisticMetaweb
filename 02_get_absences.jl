include("A0_required.jl")

# Option to use job array
# JOBARRAY = true

# Option to run for CAN
# CAN = true
if (@isdefined CAN) && CAN == true
    res = 2.5;
    occ_path = joinpath("data", "occurrences");
    pa_path = joinpath("data", "presence_absence");
    @info "Running for Canada at 2.5 arcmin resolution"
else
    res = 10.0;
    occ_path = joinpath("xtras", "occurrences");
    pa_path = joinpath("xtras", "presence_absence");
    @info "Running for Quebec at 10 arcmin resolution"
end

reference_layer = SimpleSDMPredictor(
    WorldClim, BioClim, 1; resolution = res, left=-180.0, right=-40.0, bottom=18.0, top=89.0
)

occfiles = readdir(occ_path; join=true)
filter!(!contains("all_occurrences"), occfiles) # remove backup files for QC data

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
    # jobfiles = occfiles
    trouble_species = [
        "Synaptomys_cooperi",
        "Sorex_haydeni",
        "Aplodontia_rufa",
        "Sorex_palustris",
        "Zapus_trinotatus",
        "Spilogale_putorius"
    ]
    jobfiles = string.(joinpath.(occ_path, trouble_species), ".csv")
end

ispath(pa_path) || mkpath(pa_path)

p = Progress(length(occfiles))

Threads.@threads for i in 1:length(jobfiles)
    @time try # 7- 12 min (for Canada)
        pres = similar(reference_layer, Bool)
        df = DataFrame(CSV.File(jobfiles[i]; stringtype=String))
        for r in eachrow(df)
            if !isnothing(pres[r.longitude, r.latitude])
                pres[r.longitude, r.latitude] = true
            end
        end
        Random.seed!(i)
        abs = rand(WithinRadius, pres)
        outfile = replace(replace(jobfiles[i], "occurrences" => "presence_absence"), ".csv" => ".tif")
        geotiff(outfile, [convert(Float32, replace(pres, false => nothing)), convert(Float32, replace(abs, false => nothing))])
        # next!(p)
    catch
        @warn "Error with $(jobfiles[i])"
    end
end
