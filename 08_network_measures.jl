#### Network measures ####

include("05_assemble_networks.jl")

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