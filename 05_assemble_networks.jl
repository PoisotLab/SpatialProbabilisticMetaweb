## Probabilistic distributions

# CAN = true
include("04_aggregate_sdms.jl");

D # Truncated Normal distribution per pixel

## Metaweb

# Load the previous non-reconciled metaweb if dealing with QC data
if (@isdefined CAN) && CAN == true
    mw_path = joinpath("data", "input", "canadian_thresholded_reconciled.csv")
    fit_path = joinpath("data", "input", "sdm_fit_results.csv")
else
    mw_path = joinpath("data", "input", "canadian_thresholded.csv")
    fit_path = joinpath("xtras", "input", "sdm_fit_results.csv")
end

# Parse the metaweb
mw_output = DataFrame(CSV.File(mw_path; stringtype=String))
sp = collect(keys(μ))
S = fill(0.0, length(sp), length(sp))
P = UnipartiteProbabilisticNetwork(S, sp)
for r in eachrow(mw_output)
    from = findfirst(isequal(r.from), sp)
    to = findfirst(isequal(r.to), sp)
    P[from, to] = r.score
end

# Prepare cutoff values for all species
sdm_results = CSV.read(fit_path, DataFrame)
sdm_results.species = replace.(sdm_results.species, "_" => " ")
cutoffs = Dict{String, Float64}()
for r in eachrow(sdm_results)
    cutoffs[r.species] = r.cutoff
end
# Some species have no cutoffs, so let's add an impossible one to make everything work
missing_sp = setdiff(species(P), sdm_results.species)
for m in missing_sp
    cutoffs[m] = 1.0
end
cutoffs

## Networks in space

# Get a preliminary map
function assemble_networks(
    reference_layer::SimpleSDMLayer,
    P::UnipartiteProbabilisticNetwork,
    D::Dict{String, SimpleSDMResponse},
    cutoffs::Dict{String, Float64};
    type::String="avg",
    n_itr::Int64=10,
    seed::Int64=42,
    n_regions::Tuple{Int64, Int64}=(3,3)
)
    type in ["avg", "avg_thr", "rnd", "rnd_thr"] ||
        throw(ArgumentError("type must be avg, avg_thr, rnd, or rnd_thr"))

    # Adjacency matrix
    A = adjacency(P)

    # Global BitArray to collect all values
    n_sites_total = length(reference_layer)
    networks_bit = BitArray(undef, n_sites_total, size(P)..., n_itr);

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
    for j in 1:n_lat, i in 1:n_lon
        # Subset distribution layers to the subregion
        new_coords = (
            left=lon_lims[i],
            right=lon_lims[i+1],
            bottom=lat_lims[j],
            top=lat_lims[j+1] - 0.0000000001 # temporary workaround to avoid issue with clipping coordinate
        )
        mini_reference_layer = clip(reference_layer; new_coords...)
        mini_D = Dict{String, SimpleSDMResponse}()
        for sp in String.(keys(D))
            mini_D[sp] = clip(D[sp]; new_coords...)
        end
        mini_D

        # Get corresponding indices in the Global BitArray
        inds_region = indexin(keys(mini_reference_layer), keys(reference_layer))

        # Networkify it
        sites = keys(mini_reference_layer)
        networks = zeros(Bool, length(sites), size(P)..., n_itr)
        p = Progress(length(sites))
        @threads for i in eachindex(sites)
            site = sites[i]
            s = [mini_D[s][site] for s in species(P)]
            c = [cutoffs[s] for s in species(P)]
            # Apply cooccurrence probability options
            if type == "avg"
                pcooc = @. mean(s) * mean(s)'
            elseif type == "avg_thr"
                pcooc = @. (mean(s) > c) * (mean(s) > c)'
            elseif type == "rnd"
                Random.seed!(seed)
                pcooc = @. rand(s) * rand(s)'
            elseif type == "rnd_thr"
                Random.seed!(seed + 1)
                pcooc = @. (rand(s) > c) * (rand(s) > c)'
            end
            # Extract local network and interactions
            for j in 1:size(networks, 4)
                Random.seed!(seed + j*i)
                prob_network = UnipartiteProbabilisticNetwork(pcooc .* A, species(P))
                networks[i, :, :, j] .= adjacency(rand(prob_network))
            end
            if !(@isdefined quiet) || quiet == false
                # Print progress bar
                next!(p)
            end
        end

        # Convert to BitArray to reduce memory used
        networks_bit[inds_region, :, :, :] = convert(BitArray, networks);

        # Empty memory
        networks = nothing
        mini_reference_layer = nothing
        mini_D = nothing
        GC.gc()
    end

    return networks_bit
end

# Assembly based on average
networks = assemble_networks(reference_layer, P, D, cutoffs); # 2.5 min

# Different assembly options
#=
networks_thr = assemble_networks(reference_layer, P, D, cutoffs; type="avg_thr"); # 30 sec.
networks_rnd = assemble_networks(reference_layer, P, D, cutoffs; type="rnd"); # 2 min
networks_rnd_thr = assemble_networks(reference_layer, P, D, cutoffs; type="rnd_thr"); # 30 sec.
=#
