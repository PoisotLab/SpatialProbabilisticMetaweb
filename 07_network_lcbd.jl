#### Network LCBD analysis ####

include("05_assemble_networks.jl")

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