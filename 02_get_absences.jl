include("A0_required.jl")

reference_layer = SimpleSDMPredictor(
    WorldClim, BioClim, 1; resolution = 2.5, left=-180.0, right=-40.0, bottom=18.0, top=89.0
)

occfiles = readdir(joinpath("data", "occurrences"); join=true)
filter!(contains(".csv"), occfiles) # remove subfolder

ispath(joinpath("data", "presence_absence")) || mkpath(joinpath("data", "presence_absence"))

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
    abs = rand(SurfaceRangeEnvelope, pres)
    outfile = replace(replace(occfiles[i], "occurrences" => "presence_absence"), ".csv" => ".tif")
    geotiff(outfile, [convert(Float32, replace(pres, false => nothing)), convert(Float32, replace(abs, false => nothing))])
    next!(p)
end