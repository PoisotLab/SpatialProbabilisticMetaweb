include("A0_required.jl")

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

ispath(pa_path) || mkpath(pa_path)

p = Progress(length(occfiles))

@time Threads.@threads for i in 1:length(occfiles)
    pres = similar(reference_layer, Bool)
    df = DataFrame(CSV.File(occfiles[i]; stringtype=String))
    for r in eachrow(df)
        if !isnothing(pres[r.longitude, r.latitude])
            pres[r.longitude, r.latitude] = true
        end
    end
    Random.seed!(i)
    abs = rand(WithinRadius, pres)
    outfile = replace(replace(occfiles[i], "occurrences" => "presence_absence"), ".csv" => ".tif")
    geotiff(outfile, [convert(Float32, replace(pres, false => nothing)), convert(Float32, replace(abs, false => nothing))])
    next!(p)
end