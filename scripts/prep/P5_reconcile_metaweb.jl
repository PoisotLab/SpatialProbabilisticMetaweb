#### Reconcile the species list from the metaweb with the updated GBIF taxonomy

include("../../A0_required.jl")

# Download and parse the metaweb
mw_path = joinpath("data", "input", "canadian_thresholded.csv")
mw_output = DataFrame(CSV.File(mw_path; stringtype=String))

# Turn the metaweb into a network
sp = GBIF.taxon.(unique(vcat(mw_output.from, mw_output.to)))

# Assemble GBIF and metaweb name
spdf = DataFrame(sp)
spdf.speciesid = [sp.first for sp in spdf.species]
select!(spdf, [:name, :speciesid])

# Assemble as dictionary
spdict = Dict()
for i in 1:nrow(spdf)
    spdict[spdf.name[i]] = spdf.speciesid[i]
end
spdict

# Update metaweb species with GBIF names
mw_output.from = [spdict[from] for from in mw_output.from]
mw_output.to = [spdict[to] for to in mw_output.to]

# Find duplicated rows
@rtransform!(mw_output, :pair = Pair(:from, :to))
mw_output.duplicated = [sum(p .== mw_output.pair) > 1 for p in mw_output.pair]

# Select the highest interaction score only
mw_updated = @chain begin
    # Select highest score
    @rsubset(mw_output, :duplicated == true)
    sort([:from, :to])
    groupby(:pair)
    @combine(:score = maximum(:score))
    @rtransform(:from = :pair.first, :to=:pair.second, :duplicated=true)
    select(names(mw_output))
    # Combine with non duplicated
    vcat(@rsubset(mw_output, :duplicated == false))
    select(Not([:pair, :duplicated]))
end

# Export updated metaweb
CSV.write(joinpath("data", "input", "canadian_thresholded_reconciled.csv"), mw_updated)