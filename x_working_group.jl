#### Generate data for the working group ####

include("A0_required.jl")

## Check what the data looks like
# @load joinpath("xtras", "jld2", "network_layers.jld2") layer layer_thr layer_rnd layer_rnd_thr # 350 MB file
# @load joinpath("xtras", "jld2", "network_simulations.jld2") networks networks_thr networks_rnd networks_rnd_thr # 60 MB file

## Generate the Spatial Probabilistic Metaweb
# As a layer in a JLD2 file
spatial_metaweb = layer
@save joinpath("xtras", "jld2", "working_group", "spatial_probabilistic_metaweb.jld2") spatial_metaweb

# As a 3D matrix in a JLD2 file
# Which we then convert as a layer

## Generate the co-occurrence data
# Load layers with occurrence probabilities
include("04_aggregate_sdms.jl")

# Combine in single stack
occurrence_probabilities = map(l -> broadcast(mean, l), collect(values(D)))
geotiff(joinpath("xtras", "jld2", "working_group", "occurrence_probabilities.tif"), occurrence_probabilities)

# Matrix with co-occurrence probabilities
sites = keys(reference_layer)
_nsites = length(sites)
_nspecies = length(occurrence_probabilities)
cooccurrence = zeros(Float64, _nspecies, _nspecies, _nsites)
for i in eachindex(sites)
    site = sites[i]
    sp = [occ_prob[site] for occ_prob in occurrence_probabilities]
    pcooc = sp * sp'
    cooccurrence[:, :, i] = pcooc
end
cooccurrence

# Cooccurrence as a vector of sites
cooccurrence_stack = [cooccurrence[:, :, i] for i in 1:_nsites]

# Stack of layers with co-occurrence probabilities for every pairwise combination
Base.zero(::Type{Matrix{T}}) where T = Matrix(zeros(T, (2,2)))
cooccurrence_layer = similar(reference_layer, Matrix{Float64})
cooccurrence_layer[keys(reference_layer)] = cooccurrence_stack

# Pairwise combinations
cooccurrence_vec = reduce(hcat, vec.(cooccurrence_stack))

# Export objects (2.5 GB each)
# @save joinpath("xtras", "jld2", "working_group", "cooccurrence_3D.jld2") cooccurrence
# @save joinpath("xtras", "jld2", "working_group", "cooccurrence_layer.jld2") cooccurrence_layer
# @save joinpath("xtras", "jld2", "working_group", "cooccurrence_vec.jld2") cooccurrence_vec

## Generate the uncertainty layer
occurrence_uncertainties = map(l -> broadcast(var, l), collect(values(D)))
geotiff(joinpath("xtras", "jld2", "working_group", "occurrence_uncertainties.tif"), occurrence_uncertainties)

## Export species list
CSV.write(joinpath("xtras", "jld2", "working_group", "species_list.csv"), DataFrame(sp = sp), writeheader=false)