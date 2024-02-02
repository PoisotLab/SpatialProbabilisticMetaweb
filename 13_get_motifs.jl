#### Network motifs ####

# CAN = true
# JOBARRAY = true
# MOTIF = :S4

# Load the corresponding results if dealing with CAN data or minimal example
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

# Load local networks
if (@isdefined networks) && networks isa BitArray{4}
    @info "Object networks is already defined. Not re-running previous script."
else
    @info "Running 05_assemble_networks.jl"
    include("05_assemble_networks.jl");
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
layer = SimpleSDMResponse(_mat, reference_layer)

## Motifs

# Define motifs to evaluate
if !(@isdefined MOTIF)
    MOTIF = :S4
end
SX = getproperty(unipartitemotifs(), MOTIF)

# Create a mini layer for New Brunswick if running minimal example
# CAN = true
if !(@isdefined CAN) || CAN == false
    nb_ref = read_geotiff(joinpath("data", "input", "newbrunswick_ref_10.tif"), SimpleSDMPredictor)
    nb_layer = clip(layer; boundingbox(nb_ref)...)
    _non_nb_keys = setdiff(keys(nb_layer), keys(nb_ref))
    nb_layer[_non_nb_keys] = fill(nothing, length(_non_nb_keys))
    layer = nb_layer
    GC.gc()
end

# Create layers for the motifs
SX_layer = similar(layer, Float64)
p = Progress(length(layer), "Computing motif $MOTIF")
@threads for k in keys(layer)
    SX_layer[k] = first(expected_motif_count(find_motif(layer[k], SX)))
    if !(@isdefined quiet) || quiet == false
        # Print progress bar
        next!(p)
    end
end

# Export the motifs layers
motifs_path = joinpath(results_path, "motifs")
isdir(motifs_path) || mkpath(motifs_path)
if (@isdefined JOBARRAY) && JOBARRAY == true
    write_geotiff(joinpath(motifs_path, "$MOTIF-$(string(_jobid; pad=4)).tif"), SX_layer)
else
    write_geotiff(joinpath(motifs_path, "$MOTIF.tif"), SX_layer)
end
