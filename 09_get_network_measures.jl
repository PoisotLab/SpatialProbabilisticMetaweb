#### Network measures ####

# CAN = true

# Load the corresponding sdm results if dealing with QC or CAN data
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

## Species proportions

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

# Apply on complete layer
sites = keys(layer)
T = similar(layer, Float64)
I = similar(layer, Float64)
B = similar(layer, Float64)
p = Progress(length(sites), "Species proportions")
@threads for s in sites
    T[s], I[s], B[s] = species_proportions(layer[s])
    if !(@isdefined quiet) || quiet == false
        next!(p)
    end
end

# Export layers
write_geotiff(joinpath(results_path, "proportions_T.tif"), T)
write_geotiff(joinpath(results_path, "proportions_I.tif"), I)
write_geotiff(joinpath(results_path, "proportions_B.tif"), B)
