#### Network LCBD analysis ####

include("A0_required.jl")

# Define reference layer
spatialrange = (left=-80.0, right=-50.0, bottom=45.0, top=65.0)
reference_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; spatialrange...)

# Load networks
@load joinpath("data", "jld2", "network_simulations.jld2") networks networks_thr networks_rnd networks_rnd_thr

# Get the networks LCBD
function networks_to_layer(networks::Array{Bool, 4}, reference_layer::SimpleSDMLayer)

    # Get non-zero interactions
    valued_interactions = findall(!iszero, sum(networks; dims=(1,4))[1,:,:])
    # Sum over all iterations
    by_site = sum(networks; dims=(4))

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
    lcbd_networks[keys(lcbd_networks)] = LCBD(hellinger(Zfull))[1]

    return lcbd_networks
end
lcbd_networks = networks_to_layer(networks, reference_layer)
lcbd_networks_thr = networks_to_layer(networks_thr, reference_layer)
lcbd_networks_rnd = networks_to_layer(networks_rnd, reference_layer)
lcbd_networks_rnd_thr = networks_to_layer(networks_rnd_thr, reference_layer)

# Export
geotiff(joinpath("data", "results", "lcbd_networks_mean.tif"), lcbd_networks)
geotiff(joinpath("data", "results", "lcbd_networks_rand.tif"), lcbd_networks_rnd)
geotiff(joinpath("data", "results", "lcbd_networks_mean_thr.tif"), lcbd_networks_thr)
geotiff(joinpath("data", "results", "lcbd_networks_rand_thr.tif"), lcbd_networks_rnd_thr)
