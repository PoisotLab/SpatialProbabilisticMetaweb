## Probabilistic distributions

include("04_aggregate_sdms.jl")

D # Truncated Normal distribution per pixel

## Metaweb

# Download and parse the metaweb
mw_path = joinpath("data", "input", "canadian_thresholded.csv")
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
A = zeros(Float64, richness(P), richness(P))
for (i,si) in enumerate(species(P))
    for (j,sj) in enumerate(species(P))
        A[i,j] = P[si,sj]
    end
end

# Prepare cutoff values for all species
sdm_results = CSV.read(joinpath("data", "input", "sdm_fit_results.csv"), DataFrame)
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
# networks = assemble_networks(reference_layer, P, D, A, cutoffs); # 2 min

# Different assembly options
networks_thr = assemble_networks(reference_layer, P, D, A, cutoffs; type="avg_thr"); # 30 sec.
# networks_rnd = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd"); # 2 min
# networks_rnd_thr = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd_thr"); # 30 sec.

## Network layer

# Work on the networks_thr object for now
networks = networks_thr

# Create empty objects
_empty_mat = zeros(Float64, size(networks)[2:3])
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

# Transform into layer
_mat = fill(nothing, size(reference_layer.grid))
_mat = convert(Matrix{Union{Nothing, UnipartiteProbabilisticNetwork{Float64}}}, _mat)
_inds = findall(!isnothing, reference_layer.grid)
_mat[_inds] = networks_vec
# layer = SimpleSDMResponse(_mat) # crashes, StackOverflowError
# And VS Code terminal hangs forever on StackOverflowError 😠
# Julia in a regular terminal doesn't though

## What's going on?

# Testing in a regular terminal session to avoid hanging on StackOverflowError

# Test on smaller objects
# SimpleSDMResponse(_mat[1:2, 1:2]) # still crashes
SimpleSDMResponse([fill(_mat[1], 2, 2); nothing nothing])
SimpleSDMResponse([fill(_mat[2], 2, 2); nothing nothing])
SimpleSDMResponse([_mat[1,1] _mat[1,2]; _mat[2,1] _mat[2,2]; nothing nothing])
# But those work...

# Does converting work?
convert(Matrix{Union{Nothing,eltype(_mat)}}, _mat) # meh

# Maybe defining zero type will help?
Base.zero(::Type{UnipartiteProbabilisticNetwork{T}}) where T = UnipartiteProbabilisticNetwork(zeros(T, (2,2)))
zero(UnipartiteProbabilisticNetwork{Float64})
Base.zero(::Type{UnipartiteProbabilisticNetwork{T, String}}) where T = zero(UnipartiteProbabilisticNetwork{T})
zero(UnipartiteProbabilisticNetwork{Float64, String})

# Attempt with similar
# similar(reference_layer, UnipartiteProbabilisticNetwork{Float64}) # nope
@which similar(reference_layer, UnipartiteProbabilisticNetwork{Float64})

# Decomposing similar
layer = reference_layer
TC = eltype(_mat)
TC = UnipartiteProbabilisticNetwork{Float64}
zeros(TC, size(reference_layer))
emptygrid = convert(Matrix{Union{Nothing,TC}}, zeros(TC, size(layer)))
emptygrid[findall(isnothing, layer.grid)] .= nothing
@which SimpleSDMResponse(emptygrid, layer.left, layer.right, layer.bottom, layer.top) # hangs here

# Decomposing SimpleSDMResponse
grid = emptygrid
T = eltype(grid)
l, r, b, t = layer.left, layer.right, layer.bottom, layer.top
convert(Matrix{Union{Nothing,T}}, grid)
@which SimpleSDMResponse(convert(Matrix{Union{Nothing,T}}, grid), l, r, b, t) # hangs here

# Re-attempt on smaller object
minimat = _mat[1:2, 1:2]
minimat = convert(Matrix{UnipartiteProbabilisticNetwork{Float64}}, minimat)
SimpleSDMResponse(minimat) nope

# Let's check the exact methods
_tmp_mat1 = [_mat[1,1] _mat[1,2]; _mat[2,1] _mat[2,2]; nothing nothing]
_tmp_mat2 = _mat[1:2, 1:2]
@which SimpleSDMResponse(_tmp_mat1) # works, function at line 60
@which SimpleSDMResponse(_tmp_mat2) # doesn't work, function at line 70
# why??

# Let's check the types then
typeof(_tmp_mat1)
typeof(_tmp_mat2)
eltype(_tmp_mat1)
eltype(_tmp_mat2)
# So the difference is just that _tmp_mat1 has String??

# Let's try to fix it
_tmp_mat3 = convert(Matrix{Union{Nothing, UnipartiteProbabilisticNetwork{Float64, String}}}, _mat)[1:2, 1:2]
typeof(_tmp_mat3) == typeof(_tmp_mat1)
@which SimpleSDMResponse(_tmp_mat3)
SimpleSDMResponse(_tmp_mat3)
# 🤦 of course that was just it

# Let's try it on the full matrix then
_mat2 = convert(Matrix{Union{Nothing, UnipartiteProbabilisticNetwork{Float64, String}}}, _mat)
@which SimpleSDMResponse(_mat2)
SimpleSDMResponse(_mat2)
layer = SimpleSDMResponse(_mat2, reference_layer)
conn_layer = broadcast(connectance, layer)
plot(conn_layer; c=:cividis)
# IT WORKS AAAAAAAAAAAAHHHH !!!

## LCBD values

# Get the networks LCBD
function networks_to_layer(
    networks::Array{Bool, 4}, reference_layer::SimpleSDMLayer; type::String="lcbd"
)
    type in ["lcbd", "links", "std", "mean"] ||
        throw(ArgumentError("type must be lcbd, links or probability"))

    # Get non-zero interactions
    valued_interactions = findall(!iszero, sum(networks; dims=(1,4))[1,:,:])
    # Sum over all iterations
    if type == "std"
        by_site = std(networks; dims=4)
    elseif type == "mean"
        by_site = mean(networks; dims=4)
    else
        by_site = sum(networks; dims=(4))
    end

    # Create a site x non-zero interactions matrix
    Z = zeros(eltype(by_site), (length(reference_layer), length(valued_interactions)))
    # Number of iterations where the interaction is realized
    sites = keys(reference_layer)
    for i in 1:length(sites)
        Z[i,:] = by_site[i,valued_interactions]
    end
    # Remove sites without interactions
    Zfull = Z[findall(!iszero, vec(sum(Z; dims=2))), :]

    # Network LCBD
    lcbd_networks = similar(reference_layer)
    lcbd_networks[keys(reference_layer)] = sum(Z; dims=2)
    replace!(lcbd_networks, 0 => nothing)
    if type == "lcbd"
        lcbd_networks[keys(lcbd_networks)] = LCBD(hellinger(Zfull))[1]
    elseif type == "links"
        lcbd_networks[keys(lcbd_networks)] = vec(sum(Z .> 0; dims=2)) ./ size(Z, 2)
    end

    return lcbd_networks
end
lcbd_networks = networks_to_layer(networks, reference_layer)
lcbd_networks_thr = networks_to_layer(networks_thr, reference_layer)
lcbd_networks_rnd = networks_to_layer(networks_rnd, reference_layer)
lcbd_networks_rnd_thr = networks_to_layer(networks_rnd_thr, reference_layer)
L = networks_to_layer(networks, reference_layer; type="links")
L_mean = networks_to_layer(networks, reference_layer, type="mean")
L_std = networks_to_layer(networks, reference_layer, type="std")

# Export
geotiff(joinpath("data", "results", "lcbd_networks_mean.tif"), lcbd_networks)
geotiff(joinpath("data", "results", "lcbd_networks_rand.tif"), lcbd_networks_rnd)
geotiff(joinpath("data", "results", "lcbd_networks_mean_thr.tif"), lcbd_networks_thr)
geotiff(joinpath("data", "results", "lcbd_networks_rand_thr.tif"), lcbd_networks_rnd_thr)
geotiff(joinpath("data", "results", "links.tif"), L)
geotiff(joinpath("data", "results", "links_mean.tif"), L_mean)
geotiff(joinpath("data", "results", "links_std.tif"), L_std)