#### Network LCBD analysis ####

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
    include("05_assemble_networks.jl"); # load networks
end

# Verify loaded objects
size(networks) # local networks
P # probabilistic metaweb
reference_layer # reference_layer

# Get the networks LCBD
function networks_to_layer(networks::BitArray{4}, P, reference_layer::SimpleSDMLayer)

    # Get non-zero interactions
    inds_possible = findall(!iszero, adjacency(P))
    # Sum interactions over all iterations
    freq_by_site = dropdims(mapreduce(Int8, +, networks, dims=4); dims=4);

    # Create a site x non-zero interactions matrix
    Z = @view freq_by_site[:, inds_possible]
    # Remove species without interactions
    inds_full_sp = findall(!iszero, vec(sum(Z; dims=1)))
    # Remove sites without interactions
    inds_full_sites = findall(!iszero, vec(sum(Z; dims=2)))
    # Smaller matrix
    Zfull = @view Z[inds_full_sites, inds_full_sp]

    # Network LCBD
    lcbd_networks = similar(reference_layer)
    lcbd_networks[keys(reference_layer)] = sum(Z; dims=2)
    replace!(lcbd_networks, 0 => nothing)
    lcbd_networks[keys(lcbd_networks)] = LCBD(hellinger(Zfull))[1]

    return lcbd_networks
end
lcbd_networks = networks_to_layer(networks, P, reference_layer)
# lcbd_networks_thr = networks_to_layer(networks_thr, P, reference_layer)
# lcbd_networks_rnd = networks_to_layer(networks_rnd, P, reference_layer)
# lcbd_networks_rnd_thr = networks_to_layer(networks_rnd_thr, P, reference_layer)

# Export
isdir(results_path) || mkpath(results_path)
write_geotiff(joinpath(results_path, "lcbd_networks_mean.tif"), lcbd_networks)
# geotiff(joinpath(results_path, "lcbd_networks_rand.tif"), lcbd_networks_rnd)
# geotiff(joinpath(results_path, "lcbd_networks_mean_thr.tif"), lcbd_networks_thr)
# geotiff(joinpath(results_path, "lcbd_networks_rand_thr.tif"), lcbd_networks_rnd_thr)
