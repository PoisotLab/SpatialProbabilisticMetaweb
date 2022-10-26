include("A0_required.jl")

reference_layer = SimpleSDMPredictor(
    WorldClim, BioClim, 1; left=-180.0, right=-40.0, bottom=18.0, top=89.0
)

occfiles = readdir(joinpath("data", "occurrences", "subset"); join=true)
filter!(!contains("all_occurrences"), occfiles)

ispath(joinpath("data", "presence_absence")) || mkpath(joinpath("data", "presence_absence"))

# p = Progress(length(occfiles))

# @time Threads.@threads for i in 1:length(occfiles)
@time for i in 1:length(occfiles)
    @info "Running $(occfiles[i])"
    @time begin
    pres = similar(reference_layer, Bool)
    df = DataFrame(CSV.File(occfiles[i]; stringtype=String))
    for r in eachrow(df)
        if !isnothing(pres[r.longitude, r.latitude])
            pres[r.longitude, r.latitude] = true
        end
    end
    abs = rand(SurfaceRangeEnvelope, pres)
    outfile = replace(replace(occfiles[i], "occurrences/subset" => "presence_absence"), ".csv" => ".tif")
    geotiff(outfile, [convert(Float32, replace(pres, false => nothing)), convert(Float32, replace(abs, false => nothing))])
    # next!(p)
    # GC.gc() # Won't hurt I guess? Does slow down a lot, like 54 sec. vs 14 sec.
    end
end

function get_absences(i::Int64)
    @info "Running $(occfiles[i])"
    @time begin
    pres = similar(reference_layer, Bool)
    df = DataFrame(CSV.File(occfiles[i]; stringtype=String))
    for r in eachrow(df)
        if !isnothing(pres[r.longitude, r.latitude])
            pres[r.longitude, r.latitude] = true
        end
    end
    abs = rand(SurfaceRangeEnvelope, pres)
    outfile = replace(replace(occfiles[i], "occurrences/subset" => "presence_absence"), ".csv" => ".tif")
    geotiff(outfile, [convert(Float32, replace(pres, false => nothing)), convert(Float32, replace(abs, false => nothing))])
    # next!(p)
    # GC.gc() # Won't hurt I guess? Does slow down a lot, like 54 sec. vs 14 sec.
    end
end
# @time get_absences(1)
# @time begin
# @time precompile(get_absences, ((Int64,)))
@time Threads.@threads for i in 1:length(occfiles)
# @time for i in 1:length(occfiles)
    get_absences(i)
end
# end

# - Gc.gc() slows down a lot (without threads), like 54 sec. vs 14 sec.
# And it just doesn't work with threads
# - next!(p) slows down a little, ~ 1 sec. per species
# - DataFrame(CSV.File) is the longest to precompile despite the sysimage
# - CSV.read() doesn't improve
# - stringtype=String doesn't change it
# ----
# - On a single-thread loop, the first iteration compiles everything (10s), then
# other iterations are already compiled (14s for remaining -- 24s in total)
# - Multi-thread (precompiled) gets it down to 8s
# - On a multi-thread loop, all first iterations compile everything (10s each),
# then other iterations are already compile (14s for remaining -- 23s in total)
# - So not much gain with threads on single call because on compiling
# ----
# - Custom function get_absences: 0.12s to precompile, 0.10s after, ~15s in total
# And 8s for very first compiling, 23s in total for first call
# - So no gain with multi-thread **in this case** since we run once and
# compiling time >> run time (which is very short)
# - Would be different if the runtime was much longer
# - Probably different too if the ration compile:run was different but not sure
# ----
# - @time precompile(get_absences, ((Int64,))) returns true in 0.00004s
# and doesn't precompile or speedup the loop
# - No difference with threads