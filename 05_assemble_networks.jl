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
networks = assemble_networks(reference_layer, P, D, A, cutoffs); # 2 min

# Different assembly options
# networks_thr = assemble_networks(reference_layer, P, D, A, cutoffs; type="avg_thr"); # 30 sec.
# networks_rnd = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd"); # 2 min
# networks_rnd_thr = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd_thr"); # 30 sec.

## Network layer

# Work on the networks_thr object for now
# networks = networks_thr

# Define some zero types
Base.zero(::Type{UnipartiteProbabilisticNetwork{T}}) where T = UnipartiteProbabilisticNetwork(zeros(T, (2,2)))
Base.zero(::Type{UnipartiteProbabilisticNetwork{T, String}}) where T = zero(UnipartiteProbabilisticNetwork{T})

# Define function
function network_layer(networks)
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

    # Transform into layer (option 1)
    # layer = similar(reference_layer, UnipartiteProbabilisticNetwork{Float64, String})
    layer = similar(reference_layer, eltype(networks_vec))
    layer[keys(layer)] = networks_vec

    # Transform into layer (option 2)
    _mat = fill(nothing, size(reference_layer.grid));
    _mat = convert(Matrix{Union{Nothing, eltype(networks_vec)}}, _mat);
    _inds = findall(!isnothing, reference_layer.grid);
    _mat[_inds] = networks_vec;
    layer = SimpleSDMResponse(_mat, reference_layer);

    return layer
end

# Convert all options
layer = network_layer(networks)
# layer_thr = network_layer(networks_thr)
# layer_rnd = network_layer(networks_rnd)
# layer_rnd_thr = network_layer(networks_rnd_thr)
# layers_all = [layer, layer_thr, layer_rnd, layer_rnd_thr]

# Get some network measures
Co = broadcast(connectance, layer)
# o = broadcast(omnivory, layer)
# tl = broadcast(trophic_level, layer)
L = broadcast(links, layer)
Lv = broadcast(links_var, layer)
Ld = broadcast(linkage_density, layer)
# m = broadcast(find_motif, layer)

# Plot 'em
plot(
    plot(Co; c=:cividis, title="Connectance"),
    plot(L; c=:cividis, title="Links"),
    plot(Lv; c=:cividis, title="Link variance"),
    plot(Ld; c=:cividis, title="Linkage density"),
)
savefig("figures/network_things.png")

# Checking network entropy
H = broadcast(EcologicalNetworks.entropy, layer)
plot(H; c=:cividis, title="Entropy")

# Information decomposition
H_D = broadcast(diff_entropy_uniform, layer)
H_I = broadcast(mutual_information, layer)
H_V = broadcast(variation_information, layer)
plot(
    plot(H; c=:cividis, title="Entropy (H)"),
    plot(H_D; c=:cividis, title="Difference in entropy vs uniform (D)"),
    plot(H_I; c=:cividis, title="Mutual information (I)"),
    plot(H_V; c=:cividis, title="Variation of information (V)");
    size=(900,600)
)
savefig(joinpath("figures", "new", "links_entropy.png"))

# Trivariate visualisation
begin
    tri1 = trivariate(
        H_D,
        H_I,
        H_V;
        title="Information decomposition",
        frame=:grid,
        # quantiles=false,
        # simplex=true,
    )
    xaxis!(tri1, "Longitude")
    yaxis!(tri1, "Latitude")
    tri2 = trivariatelegend!(
        H_D,
        H_I,
        H_V;
        inset=(1, bbox(0.02, 0.05, 0.35, 0.35, :top, :right)),
        subplot=2,
        red="Difference in entropy (D)",
        green="Mutual information (I)",
        blue="Variation of information (V)",
        annotationfontsize=5,
        # quantiles=false,
        # simplex=true,
    )
end
savefig(joinpath("figures", "new", "links_entropy_trivariate.png"))

## Species proportions

# Initial tests
N = layer[1]
k = degree(N)
sum(values(k)) # equal to 2L
mean(values(k))
histogram(collect(values(k)), bins=50)
isdegenerate(N) # some species have no interactions # are they even there??
simplify(N)
degree(N, dims=2)

# Extract proportions of top, bottom and intermediate species
function species_proportions(N)
    # Compute degrees
    degrees =  DataFrame(
        species = species(N)
    )
    @rtransform!(degrees, :degree = degree(N)[:species])
    @rtransform!(degrees, :out_degree = degree(N; dims=1)[:species])
    @rtransform!(degrees, :in_degree = degree(N; dims=2)[:species])
    @rtransform!(degrees, :absent = iszero(:degree))

    # Classify species as top, bottom and intermediate
    @chain degrees begin
        @rtransform!(:top = !iszero(:out_degree) && iszero(:in_degree) && iszero(:absent))
        @rtransform!(:interm = !iszero(:out_degree) && !iszero(:in_degree) && iszero(:absent))
        @rtransform!(:basal = iszero(:out_degree) && !iszero(:in_degree) && iszero(:absent))
    end

    # Extract proportions
    T = mean(degrees.top)
    I = mean(degrees.interm)
    B = mean(degrees.basal)

    return (T = T, I = I, B = B)
end

# Single network test
@time species_proportions(N)
# @time T = broadcast(v -> species_proportions(v), layer)

# Apply on complete layer
sites = keys(layer)
p = Progress(length(sites))
T = similar(layer, Float64)
I = similar(layer, Float64)
B = similar(layer, Float64)
@threads for s in sites
    T[s], I[s], B[s] = species_proportions(layer[s])
    next!(p)
end

# Visualise results
plot(T; c=:cividis)
plot(I; c=:cividis)
plot(B; c=:cividis)

# Trivariate visualisation
begin
    tri1 = trivariate(
        T,
        I,
        B;
        title="Species proportions",
        frame=:grid,
        # quantiles=false,
        # simplex=true,
    )
    xaxis!(tri1, "Longitude")
    yaxis!(tri1, "Latitude")
    tri2 = trivariatelegend!(
        T,
        I,
        B;
        inset=(1, bbox(0.02, 0.05, 0.35, 0.35, :top, :right)),
        subplot=2,
        red="Top",
        green="Intermediate",
        blue="Basal",
        annotationfontsize=7,
        # quantiles=false,
        # simplex=true,
    )
end
savefig(joinpath("figures", "new", "links_proportions_trivariate.png"))

# ## LCBD values

# # Get the networks LCBD
# function networks_to_layer(
#     networks::Array{Bool, 4}, reference_layer::SimpleSDMLayer; type::String="lcbd"
# )
#     type in ["lcbd", "links", "std", "mean"] ||
#         throw(ArgumentError("type must be lcbd, links or probability"))

#     # Get non-zero interactions
#     valued_interactions = findall(!iszero, sum(networks; dims=(1,4))[1,:,:])
#     # Sum over all iterations
#     if type == "std"
#         by_site = std(networks; dims=4)
#     elseif type == "mean"
#         by_site = mean(networks; dims=4)
#     else
#         by_site = sum(networks; dims=(4))
#     end

#     # Create a site x non-zero interactions matrix
#     Z = zeros(eltype(by_site), (length(reference_layer), length(valued_interactions)))
#     # Number of iterations where the interaction is realized
#     sites = keys(reference_layer)
#     for i in 1:length(sites)
#         Z[i,:] = by_site[i,valued_interactions]
#     end
#     # Remove sites without interactions
#     Zfull = Z[findall(!iszero, vec(sum(Z; dims=2))), :]

#     # Network LCBD
#     lcbd_networks = similar(reference_layer)
#     lcbd_networks[keys(reference_layer)] = sum(Z; dims=2)
#     replace!(lcbd_networks, 0 => nothing)
#     if type == "lcbd"
#         lcbd_networks[keys(lcbd_networks)] = LCBD(hellinger(Zfull))[1]
#     elseif type == "links"
#         lcbd_networks[keys(lcbd_networks)] = vec(sum(Z .> 0; dims=2)) ./ size(Z, 2)
#     end

#     return lcbd_networks
# end
# lcbd_networks = networks_to_layer(networks, reference_layer)
# lcbd_networks_thr = networks_to_layer(networks_thr, reference_layer)
# lcbd_networks_rnd = networks_to_layer(networks_rnd, reference_layer)
# lcbd_networks_rnd_thr = networks_to_layer(networks_rnd_thr, reference_layer)
# L = networks_to_layer(networks, reference_layer; type="links")
# L_mean = networks_to_layer(networks, reference_layer, type="mean")
# L_std = networks_to_layer(networks, reference_layer, type="std")

# # Export
# geotiff(joinpath("data", "results", "lcbd_networks_mean.tif"), lcbd_networks)
# geotiff(joinpath("data", "results", "lcbd_networks_rand.tif"), lcbd_networks_rnd)
# geotiff(joinpath("data", "results", "lcbd_networks_mean_thr.tif"), lcbd_networks_thr)
# geotiff(joinpath("data", "results", "lcbd_networks_rand_thr.tif"), lcbd_networks_rnd_thr)
# geotiff(joinpath("data", "results", "links.tif"), L)
# geotiff(joinpath("data", "results", "links_mean.tif"), L_mean)
# geotiff(joinpath("data", "results", "links_std.tif"), L_std)