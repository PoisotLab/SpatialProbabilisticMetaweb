## Probabilistic distributions

QC = true
include("04_aggregate_sdms.jl")

D # Truncated Normal distribution per pixel

## Metaweb

# Load the previous non-reconciled metaweb if dealing with QC data
if (@isdefined QC) && QC == true
    mw_path = joinpath("data", "input", "canadian_thresholded.csv")
    results_path = joinpath("xtras", "input", "sdm_fit_results.csv")
else
    mw_path = joinpath("data", "input", "canadian_thresholded_reconciled.csv")
    results_path = joinpath("data", "input", "sdm_fit_results.csv")
end

# Parse the metaweb
mw_output = DataFrame(CSV.File(mw_path; stringtype=String))
sp = collect(keys(μ))
S = fill(0.0, length(sp), length(sp))
P = UnipartiteProbabilisticNetwork(S, sp)
for r in eachrow(mw_output)
    from = findfirst(isequal(r.from), sp)
    to = findfirst(isequal(r.to), sp)
    P[from, to] = r.score
end

# Interaction matrix
A = adjacency(P)

# Prepare cutoff values for all species
sdm_results = CSV.read(results_path, DataFrame)
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

## Networks in space

# FOR I IN REGIONS

# Divide into smaller regions
spatialrange
new_coords = (
    left=(spatialrange.left + (spatialrange.right - spatialrange.left)/2),
    top=(spatialrange.top - (spatialrange.top - spatialrange.bottom)/2)
)
mini_reference_layer = clip(reference_layer; new_coords...)
mini_D = Dict{String, SimpleSDMResponse}()
for sp in String.(keys(D))
    mini_D[sp] = clip(D[sp]; new_coords...)
end
mini_D

# Get a preliminary map
function assemble_networks(
    reference_layer::SimpleSDMLayer,
    P::UnipartiteProbabilisticNetwork,
    D::Dict{String, SimpleSDMResponse},
    A::Matrix,
    cutoffs::Dict{String, Float64};
    type::String="avg",
    n_itr::Int64=10,
    seed::Int64=42
)
    type in ["avg", "avg_thr", "rnd", "rnd_thr"] ||
        throw(ArgumentError("type must be avg, avg_thr, rnd, or rnd_thr"))

    sites = keys(reference_layer)
    networks = zeros(Bool, length(sites), size(P)..., n_itr)
    p = Progress(length(sites))
    Threads.@threads for i in eachindex(sites)
        site = sites[i]
        s = [D[s][site] for s in species(P)]
        c = [cutoffs[s] for s in species(P)]
        if type == "avg"
            pcooc = @. mean(s) * mean(s)'
        elseif type == "avg_thr"
            pcooc = @. (mean(s) > c) * (mean(s) > c)'
        elseif type == "rnd"
            Random.seed!(seed)
            pcooc = @. rand(s) * rand(s)'
        elseif type == "rnd_thr"
            Random.seed!(seed + 1)
            pcooc = @. (rand(s) > c) * (rand(s) > c)'
        end
        for j in 1:size(networks, 4)
            Random.seed!(seed + j*i)
            prob_network = UnipartiteProbabilisticNetwork(pcooc .* A, species(P))
            networks[i, :, :, j] .= adjacency(rand(prob_network))
        end
        next!(p)
    end
    return networks
end

# Assembly based on average
networks = assemble_networks(mini_reference_layer, P, mini_D, A, cutoffs); # 2 min
bitnetworks = convert(BitArray, networks);
varinfo(r"networks")
Base.summarysize(networks)

Base.summarysize(zero(Bool)) # 1
Base.summarysize(zero(Int8)) # 1
Base.summarysize(zero(Float16)) # 2
Base.summarysize(BitSet(1)) # 64


# Different assembly options
#=
networks_thr = assemble_networks(reference_layer, P, D, A, cutoffs; type="avg_thr"); # 30 sec.
networks_rnd = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd"); # 2 min
networks_rnd_thr = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd_thr"); # 30 sec.
=#

## Network layer

# Work on the networks_thr object for now
# networks = networks_thr

# Define some zero types
Base.zero(::Type{UnipartiteProbabilisticNetwork{T}}) where T = UnipartiteProbabilisticNetwork(zeros(T, (2,2)))
Base.zero(::Type{UnipartiteProbabilisticNetwork{T, String}}) where T = zero(UnipartiteProbabilisticNetwork{T})

# Define function
function network_layer(networks, reference_layer)
    # Create sparse array
    M = 100;
    N = 1000;
    nz = 2000; # number of nonzeros
    I = rand(1:M, nz); # dummy I indices
    J = rand(1:N, nz); # dummy J indices
    V = randn(nz); # dummy matrix values
    _empty_mat = sparse(I, J, V, M, N)

    # Create empty objects
    _empty_network = UnipartiteProbabilisticNetwork(_empty_mat, species(P))

    # With threads
    networks_vec = fill(_empty_network, size(networks)[1])
    # networks_vec = networks_vec[1:1000]
    @threads for i in eachindex(networks_vec)
        # Extract a single site
        local_site = @view networks[i, :, :, :];

        # Extract local interaction matrix
        local_mat = dropdims(mean(local_site; dims=ndims(local_site)), dims=ndims(local_site))

        # Transform into network
        networks_vec[i] = UnipartiteProbabilisticNetwork(local_mat, species(P))
    end
    networks_vec

    # Transform into layer (option 1)
    # layer = similar(reference_layer, eltype(networks_vec))
    # layer[keys(layer)] = networks_vec

    # Transform into layer (option 2)
    _mat = fill(nothing, size(reference_layer.grid));
    _mat = convert(Matrix{Union{Nothing, eltype(networks_vec)}}, _mat);
    _inds = findall(!isnothing, reference_layer.grid);
    _mat[_inds] = networks_vec;
    layer = SimpleSDMResponse(_mat, reference_layer);

    return layer
end

# Convert all options
layer = network_layer(networks, mini_reference_layer)
#=
layer_thr = network_layer(networks_thr)
layer_rnd = network_layer(networks_rnd)
layer_rnd_thr = network_layer(networks_rnd_thr)
=#

## Export everything to JLD2

# Export
#=
isdir(joinpath("xtras", "jld2")) || mkpath(joinpath("xtras", "jld2"))
@save joinpath("xtras", "jld2", "network_layers.jld2") layer layer_thr layer_rnd layer_rnd_thr
@save joinpath("xtras", "jld2", "network_simulations.jld2") {compress=true} networks networks_thr networks_rnd networks_rnd_thr
=#