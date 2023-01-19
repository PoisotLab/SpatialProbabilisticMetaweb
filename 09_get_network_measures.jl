#### Network measures ####

# QC = true
include("A0_required.jl")

# Load the previous sdm results if dealing with QC data
if (@isdefined QC) && QC == true
    results_path = joinpath("xtras", "results")
else
    results_path = joinpath("data", "results")
end

# Load networks
include("05_assemble_networks.jl");

## Network layer

# Define some zero types
Base.zero(::Type{UnipartiteProbabilisticNetwork{T}}) where T = UnipartiteProbabilisticNetwork(zeros(T, (2,2)))
Base.zero(::Type{UnipartiteProbabilisticNetwork{T, String}}) where T = zero(UnipartiteProbabilisticNetwork{T})

# Define function
function network_layer(networks, reference_layer; FloatType::DataType=Float64)
    # Create empty objects
    _empty_mat = zeros(FloatType, size(networks)[2:3])
    _empty_network = UnipartiteProbabilisticNetwork(_empty_mat, species(P))

    # With threads
    networks_vec = fill(_empty_network, size(networks)[1])
    @threads for i in eachindex(networks_vec)
        # Extract a single site
        local_site = @view networks[i, :, :, :];

        # Extract local interaction matrix
        local_mat = dropdims(mean(local_site; dims=ndims(local_site)), dims=ndims(local_site));
        local_mat = convert(Matrix{FloatType}, local_mat);

        # Transform into network
        networks_vec[i] = UnipartiteProbabilisticNetwork(local_mat, species(P))
    end
    networks_vec

    # Transform into layer
    _mat = fill(nothing, size(reference_layer.grid));
    _mat = convert(Matrix{Union{Nothing, eltype(networks_vec)}}, _mat);
    _inds = findall(!isnothing, reference_layer.grid);
    _mat[_inds] = networks_vec;
    layer = SimpleSDMResponse(_mat, reference_layer);

    return layer
end

# Convert all options
layer = network_layer(networks, reference_layer)
# @time layer16 = network_layer(networks, reference_layer; FloatType=Float16)

## Network properties

# Get some network measures
Co = broadcast(connectance, layer)
L = broadcast(links, layer)
Lv = broadcast(links_var, layer)
Ld = broadcast(linkage_density, layer)
# o = broadcast(omnivory, layer)
# tl = broadcast(trophic_level, layer)
# m = broadcast(find_motif, layer)

# Export the link layers
geotiff(joinpath(results_path, "links_mean.tif"), L)
geotiff(joinpath(results_path, "links_var.tif"), Lv)

# Plot 'em
plot(
    plot(Co; c=:cividis, title="Connectance"),
    plot(L; c=:cividis, title="Links"),
    plot(Lv; c=:cividis, title="Link variance"),
    plot(Ld; c=:cividis, title="Linkage density");
    size=(900,600)
)