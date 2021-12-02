## What's needed

# Load aggregation script
include("04_aggregate_sdms.jl")

# Prepare cutoff values for all species
sdm_results = CSV.read("sdm_fit_results.csv", DataFrame)
sdm_results.species = replace.(sdm_results.species, "_" => " ")
cutoffs = Dict{String, Float64}()
for r in eachrow(sdm_results)
    cutoffs[r.species] = r.cutoff
end
# Some species have no cutoffs, so let's add an impossible one to make everything work
missing_sp = setdiff(species(P), sdm_results.species)
for m in missing_sp
    cutoffs[m] = 1.0
end
cutoffs

## Getting to the 4 richness options

# Rand should be similar process as mean
sp1 = "Tamias striatus"
Dmean = broadcast(mean, D[sp1])
Drand = broadcast(rand, D[sp1])

# Compare them
plot(
    plot(Dmean),
    plot(Drand),
    layout=(2,1),
    size=(600,600),
)

# Compare with μ
plot(
    plot(μ[sp1]),
    plot(Dmean),
    layout=(2,1),
    size=(600,600),
)

# Why is μ not like mean??
mean(μ[sp1])
mean(Dmean)
# ??

# Is it the number of values greater than zero?
using DataFramesMeta
@chain begin μ[sp1]
    collect
    filter(>=(0), _)
end
@chain begin Dmean
    collect
    filter(>=(0), _)
end
# Not quite

# Is it because of the truncated normal?
dnorm = Normal(μ[sp1][1], σ[sp1][1])
dtrunc = Truncated(dnorm, 0, 1)
mean(dnorm)
mean(dtrunc)
std(dnorm)
std(dtrunc)
# ah-ah!

# Plot them to see the difference
Drand_cut = convert(Float32, broadcast(>(cutoffs[sp1]), Drand))
Dmean_cut = convert(Float32, broadcast(>(cutoffs[sp1]), Dmean))
plot(
    plot(Dmean_cut),
    plot(Drand_cut),
    layout=(2,1),
    size=(600,600),
)
