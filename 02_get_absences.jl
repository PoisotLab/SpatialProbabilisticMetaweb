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
# if (@isdefined JOBARRAY) && JOBARRAY == true
    _jobid = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
    _jobcount = parse(Int64, get(ENV, "SLURM_ARRAY_TASK_COUNT", "1"))

    _ntot = length(occfiles)
    _nfiles = ceil(Int64, _ntot/_jobcount)
    _idx = [[(i*_nfiles + 1) for i in 0:(_jobcount-1)]..., _ntot]

    _idx1 = _idx[_jobid]
    _idx2 = _idx[_jobid+1]

    jobfiles = occfiles[_idx1:_idx2]

    # How inequal is the filesize between jobs?
    gp0 = repeat(1:_jobcount, inner=_nfiles)[1:_ntot]
    gp0_f = [occfiles[findall(==(i), gp0)] for i in 1:_jobcount]
    gp0_fs = [sum(filesize.(f))/1000 for f in gp0_f]
    extrema(gp0_fs)

    # Balance larger files between jobs
    gp = repeat(vcat(collect(1:_jobcount), collect(_jobcount:-1:1)), Int(_nfiles/2))[1:_ntot]
    gp_f = [sort(occfiles, by=filesize, rev=true)[findall(==(i), gp)] for i in 1:_jobcount]
    gp_fs = [sum(filesize.(f))/1000 for f in gp_f]
    extrema(gp_fs)

    # Balance files according to filesize
    fs = filesize.(occfiles)
    fs_rel = fs/sum(fs)
    fs_cumu = cumsum(sort(fs_rel, rev=true))
    fs_prop = fs_cumu/(1/_jobcount)
    gp2 = round.(Int, fs_prop)
    gp2_f = [sort(occfiles, by=filesize, rev=true)[findall(==(i), gp2)] for i in 1:_jobcount]
    gp2_fs = [sum(filesize.(f))/1000 for f in gp2_f]
    extrema(gp2_fs)

    # Fix imbalance between two final groups
    nmin = argmin(gp2_fs)
    nmax = argmax(gp2_fs)
    gp2_fs_max = filesize.(gp2_f[nmax])/1000
    nmax_min = argmax(gp2_fs_max)
    if maximum(gp2_fs) - minimum(gp2_fs) > maximum(gp2_fs_max)
        _to_change = gp2_f[nmax][nmax_min]
        gp2_f[nmin] = vcat(gp2_f[nmin], _to_change)
        gp2_f[nmax] = filter(!=(_to_change), gp2_f[nmax])
    end
    gp2_fs = [sum(filesize.(f))/1000 for f in gp2_f]
    extrema(gp2_fs)


# else
    obfiles = occfiles
# end

ispath(pa_path) || mkpath(pa_path)

p = Progress(length(occfiles))

Threads.@threads for i in 1:length(jobfiles)
    try
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
