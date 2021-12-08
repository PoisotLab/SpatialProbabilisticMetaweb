using SimpleSDMLayers
using DataFrames
import CSV
using ProgressMeter

reference_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; left=-180., right=-40., bottom=18., top=89.)

occfiles = readdir(joinpath("data", "occurrences"); join=true)

ispath(joinpath("data", "presence_absence")) || mkpath(joinpath("data", "presence_absence"))

p = Progress(length(occfiles))

Threads.@threads for i in 1:length(occfiles)
    pres = similar(reference_layer, Bool)
    df = DataFrame(CSV.File(occfiles[i]; stringtype=String))
    for r in eachrow(df)
        if !isnothing(pres[r.longitude, r.latitude])
            pres[r.longitude, r.latitude] = true
        end
    end
    abs = rand(SurfaceRangeEnvelope, pres)
    outfile = replace(replace(occfiles[i], "occurrences" => "presence_absence"), ".csv" => ".tif")
    geotiff(outfile, [convert(Float32, replace(pres, false => nothing)), convert(Float32, replace(abs, false => nothing))])
    next!(p)
    GC.gc() # Won't hurt I guess?
end
