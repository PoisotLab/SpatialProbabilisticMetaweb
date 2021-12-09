using SimpleSDMLayers
using EcologicalNetworks
using GBIF
using DataFrames
import CSV
using ProgressMeter

# Allow GBIF taxa as network nodes
EcologicalNetworks._check_species_validity(::Type{GBIFTaxon}) = nothing

# Download and parse the metaweb
mw_path = joinpath("data", "input", "canadian_thresholded.csv")
download("https://raw.githubusercontent.com/PoisotLab/MetawebTransferLearning/main/artifacts/canadian_thresholded.csv", mw_path)
mw_output = DataFrame(CSV.File(mw_path; stringtype=String))

# Turn the metaweb into a network
sp = GBIF.taxon.(unique(vcat(mw_output.from, mw_output.to)))
S = fill(0.0, length(sp), length(sp))
P = UnipartiteProbabilisticNetwork(S, sp)

for r in eachrow(mw_output)
    from = findfirst(t -> isequal(r.from)(t.name), sp)
    to = findfirst(t -> isequal(r.to)(t.name), sp)
    P[from, to] = r.score
end

# Occurrences for every taxon
p = Progress(richness(P))

ispath(joinpath("data", "occurrences")) || mkpath(joinpath("data", "occurrences"))

Threads.@threads for i in 1:richness(P)
    this_sp = species(P)[i]
    _spname = this_sp.name
    _file = joinpath("data", "occurrences", replace(_spname, " " => "_")*".csv")
    if !isfile(_file)
        this_occ = occurrences(this_sp, "hasCoordinate" => true, "limit" => 300, "decimalLatitude"=>"10,90", "decimalLongitude" => "-175,-45")
        while length(this_occ) < min(2_000, size(this_occ))
            occurrences!(this_occ)
        end
        df = DataFrame(name = String[], latitude = Float64[], longitude = Float64[])
        for o in this_occ
            push!(df, (
                _spname, o.latitude, o.longitude
            ))
        end
        CSV.write(_file, unique(df))
    end
    next!(p)
    GC.gc() # There is a memory leak I'm not too sure about - maybe the DataFrame?
end