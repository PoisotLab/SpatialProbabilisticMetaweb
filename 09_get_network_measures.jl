#### Network measures ####

# CAN = true
include("05_assemble_networks.jl"); # Load networks

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

## Network layer

# Verify loaded objects
size(networks) # local networks
reference_layer # reference_layer

# Define some zero types
Base.zero(::Type{UnipartiteProbabilisticNetwork{T}}) where T = UnipartiteProbabilisticNetwork(zeros(T, (2,2)))
Base.zero(::Type{UnipartiteProbabilisticNetwork{T, String}}) where T = zero(UnipartiteProbabilisticNetwork{T})

# Create empty objects
_empty_mat = zeros(Float64, size(networks)[2:3])
_empty_network = UnipartiteProbabilisticNetwork(_empty_mat, species(P))

# With threads
networks_vec = fill(_empty_network, size(networks)[1])
@threads for i in eachindex(networks_vec)
    # Extract a single site
    local_site = @view networks[i, :, :, :];

    # Extract local interaction matrix
    local_mat = dropdims(mean(local_site; dims=ndims(local_site)), dims=ndims(local_site));
    local_mat = convert(Matrix{eltype(_empty_mat)}, local_mat);

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

# Get the network measures
Co = connectance.(layer)
L = links.(layer)
Lv = links_var.(layer)
Ld = linkage_density.(layer)

# Export the measures layers
isdir(results_path) || mkpath(results_path)
write_geotiff(joinpath(results_path, "connectance.tif"), Co)
write_geotiff(joinpath(results_path, "links_mean.tif"), L)
write_geotiff(joinpath(results_path, "links_var.tif"), Lv)
write_geotiff(joinpath(results_path, "links_density.tif"), Ld)

## Motifs

# Define motifs to evaluate
S4 = unipartitemotifs().S4
S5 = unipartitemotifs().S5

# Create a mini layer for New Brunswick
# NB = true
if (@isdefined NB) && NB == true
    nb_ref = read_geotiff(joinpath("data", "input", "newbrunswick_ref_10.tif"), SimpleSDMPredictor)
    nb_layer = clip(layer; boundingbox(nb_ref)...)
    _non_nb_keys = setdiff(keys(nb_layer), keys(nb_ref))
    nb_layer[_non_nb_keys] = fill(nothing, length(_non_nb_keys))
    layer = nb_layer
    GC.gc()
end

# Create layers for the motifs
S4_layer = similar(layer, Float64)
S5_layer = similar(layer, Float64)
p = Progress(length(layer))
@time @threads for k in keys(layer)
    S4_layer[k] = first(expected_motif_count(find_motif(layer[k], S4)))
    S5_layer[k] = first(expected_motif_count(find_motif(layer[k], S5)))
    next!(p)
end

# Export the motifs layers
isdir(results_path) || mkpath(results_path)
write_geotiff(joinpath(results_path, "S4.tif"), S4_layer)
write_geotiff(joinpath(results_path, "S5.tif"), S5_layer)
