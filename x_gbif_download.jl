using SimpleSDMLayers
using GBIF
using DataFramesMeta
using StatsBase
import CSV

# Download and parse the metaweb
mw_path = joinpath("data", "input", "canadian_thresholded.csv")
mw_output = DataFrame(CSV.File(mw_path; stringtype=String))

# Turn the metaweb into a network
sp = GBIF.taxon.(unique(vcat(mw_output.from, mw_output.to)))

# Find potential duplicates
spdf = DataFrame(sp)
select!(spdf, Not([:kingdom, :phylum, :class, :order, :family, :genus]))
id = [sp.second for sp in spdf.species]
unique(id)
spdf.duplicate = [sum(sp .== spdf.species) > 1 for sp in spdf.species]
@chain begin
    spdf
    @rsubset(:duplicate == true)
    show(_, allrows=true)
end

# What happens when we have duplicates?
occurrences(GBIF.taxon("Canis rufus", strict=true), "hasCoordinate" => true)
occurrences(GBIF.taxon("Canis lupus", strict=true), "hasCoordinate" => true)
# Same occurrences

## Evaluate number of observations
# How does this work
sp1 = occurrences(sp[1], "hasCoordinate" => true)
length(sp1) # downloaded obs
size(sp1) # total obs

# Get occurrences for all species
occ = Vector{GBIFRecords}(undef, length(sp))
for i in eachindex(sp)
    occ[i] = occurrences(
        sp[i],
        "hasCoordinate" => true,
        "limit" => 1,
        "decimalLatitude"=>"10,90",
        "decimalLongitude" => "-175,-45"
    )
end
occ
occdf = select(spdf, :name)
occdf.nocc = [size(o) for o in occ]
sort(occdf, :nocc, rev=true) |> x -> show(x, allrows=true)
# Queries seem limited at 100_000 obs
# Some species have less than 100 obs

# Check the current file size
occdf.ncurrent = zeros(Int64, nrow(occdf))
for r in eachrow(occdf)
    _spname = replace(r.name, " " => "_")
    _df = CSV.read(joinpath("data", "occurrences", _spname*".csv"), DataFrame)
    r.ncurrent = nrow(_df)
end
sort(occdf, :nocc, rev=true) |> x -> show(x, allrows=true) # makes no sense...
# Ah, probably because of the coordinates!
# Still doesn't make much sense
# Can't only be because of recent increase in number of observations?
# Maybe taxonomic backbone update?

# How many more observations do we have
sum(occdf.nocc) # 1,780,025
sum(occdf.ncurrent) # 159,471

## Extract the keys
spdf.key = [sp.second for sp in spdf.species]
spdf
spkey = @chain begin
    leftjoin(occdf, spdf, on = :name)
    select(:name, :nocc, :ncurrent, :key)
    sort(:nocc, rev=true)
end
[spkey.key] # copy paste this in curl command

## Read the downloaded data

# Read a smaller dataset to precompile
df = DataFrame(CSV.File(joinpath("data/input/canadian_thresholded.csv")))

# Attempt to read the complete dataset (should all fail)
occ_path = joinpath("data", "occurrences", "subset", "all_occurrences.csv")
df = DataFrame(CSV.File(occ_path))
df = DataFrame(CSV.File(occ_path; stringtype=String))
df = DataFrame(CSV.File(occ_path, delim="\t", types=Dict(:gbifID => String)))
df = DataFrame(CSV.File(occ_path, delim="\t", types=Dict(:gbifID => String), quotechar='"'))
df = DataFrame(CSV.File(occ_path, select=[:species, :decimalLatitude, :decimalLongitude]))

# Attempt with fewer rows
df = DataFrame(CSV.File(occ_path, limit=100))
df = DataFrame(CSV.File(occ_path, limit=100, delim="\t"))
df = DataFrame(CSV.File(occ_path, limit=100, delim="\t", types=Dict(:gbifID => String)))
df = DataFrame(CSV.File(occ_path, limit=1_000)) # not reading the right number of rows...
df = DataFrame(CSV.File(occ_path, limit=11_800)) # when it stops working

# Attempt with dataset exported with R
df = DataFrame(CSV.File(joinpath("data", "occurrences", "subset", "all_occurrences2.csv"), delim="\t"))
# Not 100_000 rows??
# Same with R though...

# Attempt with problem dataset
df = DataFrame(CSV.File(joinpath("data", "occurrences", "subset", "all_occurrences_problems.csv"), delim="\t"))
