#### Network measures ####

QC = true
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

### Subregions
function divide_subregions(networks::BitArray{4}, reference_layer; n_regions::Tuple{Int64, Int64}=(3,3))
    # Define global range
    spatialrange = boundingbox(reference_layer)
    lat_range = spatialrange.top - spatialrange.bottom
    lon_range = spatialrange.right - spatialrange.left

    # Divide into smaller regions
    n_lat = n_regions[1]
    n_lon = n_regions[2]
    lat_step = floor(lat_range/n_lat)
    lon_step = floor(lon_range/n_lon)
    lat_lims = (spatialrange.bottom .+ collect(0:(n_lat-1)) .* lat_step..., spatialrange.top)
    lon_lims = (spatialrange.left .+ collect(0:(n_lon-1)) .* lon_step..., spatialrange.right)

    # Get the networks one subregion at the time
    region_coords = []
    region_inds = []
    for j in 1:n_lat, i in 1:n_lon
        # Subset distribution layers to the subregion
        _region_coords = (
            left=lon_lims[i],
            right=lon_lims[i+1],
            bottom=lat_lims[j],
            top=lat_lims[j+1]
        )
        push!(region_coords, _region_coords)
        mini_reference_layer = clip(reference_layer; _region_coords...)

        # Get corresponding indices in the Global BitArray
        _region_inds = indexin(keys(mini_reference_layer), keys(reference_layer))
        push!(region_inds, _region_inds)
    end

    return (region_coords = region_coords, region_inds = region_inds)
end

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

function get_network_properties(layer)
    # Get some network measures
    Co = broadcast(connectance, layer)
    L = broadcast(links, layer)
    Lv = broadcast(links_var, layer)
    Ld = broadcast(linkage_density, layer)
    # o = broadcast(omnivory, layer)
    # tl = broadcast(trophic_level, layer)
    # m = broadcast(find_motif, layer)

    return (Co = Co, L = L, Lv = Lv, Ld = Ld)
end

function get_network_properties(networks, reference_layer; kw...)
    region_coords, region_inds = divide_subregions(networks, reference_layer; kw...)
    region_reference_layers = [clip(reference_layer; rc...) for rc in region_coords]
    reg_Co = similar.(region_reference_layers)
    reg_L = similar.(region_reference_layers)
    reg_Lv = similar.(region_reference_layers)
    reg_Ld = similar.(region_reference_layers)
    for i in eachindex(region_coords)
        reg_network = @view networks[region_inds[i], :, :, :];
        reg_reference_layer = region_reference_layers[i]

        reg_layer = network_layer(reg_network, reg_reference_layer; FloatType = Float32)

        reg_Co[i], reg_L[i], reg_Lv[i], reg_Ld[i] = get_network_properties(reg_layer)

        reg_network = nothing
        reg_reference_layer = nothing
        reg_layer = nothing
        GC.gc()
    end

    Co = similar(reference_layer)
    L = similar(reference_layer)
    Lv = similar(reference_layer)
    Ld = similar(reference_layer)
    measures = (Co, L, Lv, Ld)
    measures_by_region = (reg_Co, reg_L, reg_Lv, reg_Ld)
    for i in eachindex(measures)
        for reg in measures_by_region[i]
            measures[i][keys(reg)] = collect(reg)
        end
    end

    return (Co = Co, L = L, Lv = Lv, Ld = Ld)
end

# Get network properties
Co, L, Lv, Ld = get_network_properties(networks, reference_layer)

# Export the link layers
geotiff(joinpath(results_path, "links_mean.tif"), L)
geotiff(joinpath(results_path, "links_var.tif"), Lv)
