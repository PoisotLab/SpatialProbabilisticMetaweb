#### Extra network measures ###

QC = true

# Load initial network measures
include("09_get_network_measures.jl")

## Entropy

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
savefig(joinpath("figures", "entropy.png"))

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
savefig(joinpath("figures", "entropy_trivariate.png"))

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
savefig(joinpath("figures", "proportions_trivariate.png"))