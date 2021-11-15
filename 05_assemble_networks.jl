## Probabilistic distributions

include("04_aggregate_sdms.jl")

D # Truncated Normal distribution per pixel
lcbd_species

## Metaweb

# Download and parse the metaweb
mw_output = DataFrame(CSV.File("canadian_thresholded.csv"; stringtype=String))
sp = collect(keys(μ))
S = fill(0.0, length(sp), length(sp))
P = UnipartiteProbabilisticNetwork(S, sp)
for r in eachrow(mw_output)
    from = findfirst(isequal(r.from), sp)
    to = findfirst(isequal(r.to), sp)
    P[from, to] = r.score
end

# Interaction matrix
A = zeros(Float64, richness(P), richness(P))
for (i,si) in enumerate(species(P))
    for (j,sj) in enumerate(species(P))
        A[i,j] = P[si,sj]
    end
end

# Prepare cutoff values for all species
sdm_results = CSV.read("sdm_fit_results.csv", DataFrame)
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
    A::Matrix,
    cutoffs::Dict{String, Float64};
    type::String="avg",
    n_itr::Int64=10,
    seed::Int64=42
)
    type in ["avg", "avg_thr", "rnd", "rnd_thr"] ||
        throw(ArgumentError("type must be avg, avg_thr, rnd, or rnd_thr"))

    sites = keys(reference_layer)
    networks = zeros(Bool, length(sites), size(P)..., n_itr)
    p = Progress(length(sites))
    Threads.@threads for i in eachindex(sites)
        site = sites[i]
        s = [D[s][site] for s in species(P)]
        c = [cutoffs[s] for s in species(P)]
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
        for j in 1:size(networks, 4)
            Random.seed!(seed + j*i)
            prob_network = UnipartiteProbabilisticNetwork(pcooc .* A, species(P))
            networks[i, :, :, j] .= adjacency(rand(prob_network))
        end
        next!(p)
    end
    return networks
end

# Assembly based on average
networks = assemble_networks(reference_layer, P, D, A, cutoffs); # 2 min

# Different assembly options
networks_thr = assemble_networks(reference_layer, P, D, A, cutoffs; type="avg_thr"); # 30 sec.
networks_rnd = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd"); # 2 min
networks_rnd_thr = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd_thr"); # 30 sec.

## LCBD values

# Get the networks LCBD
function networks_to_layer(
    networks::Array{Bool, 4}, reference_layer::SimpleSDMLayer; type::String="lcbd"
)
    type in ["lcbd", "links"] ||
        throw(ArgumentError("type must be lcbd or links"))

    # Get non-zero interactions
    valued_interactions = findall(!iszero, sum(networks; dims=(1,4))[1,:,:])
    # Sum over all iterations
    by_site = sum(networks; dims=(4))

    # Create a site x non-zero interactions matrix
    Z = zeros(Int64, (length(reference_layer), length(valued_interactions)))
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
    else
        lcbd_networks[keys(lcbd_networks)] = vec(sum(Z .> 0; dims=2)) ./ size(Z, 2)
    end

    return lcbd_networks
end
lcbd_networks = networks_to_layer(networks, reference_layer)
lcbd_networks_thr = networks_to_layer(networks_thr, reference_layer)
lcbd_networks_rnd = networks_to_layer(networks_rnd, reference_layer)
lcbd_networks_rnd_thr = networks_to_layer(networks_rnd_thr, reference_layer)
links = networks_to_layer(networks, reference_layer; type="links")

## Bivariate maps

# Prepare bivariate colors
p0 = colorant"#e8e8e8"
bv_pal_4 = (p0=p0, p1=colorant"#be64ac", p2=colorant"#5ac8c8")

# Bivariate LCBD
biv_plots = []
for lcbd_n in [lcbd_networks, lcbd_networks_thr, lcbd_networks_rnd, lcbd_networks_rnd_thr]
    bp = bivariate(lcbd_n, lcbd_species; quantiles=true, bv_pal_4..., classes=3)
    bp = bivariatelegend!(
        lcbd_n,
        lcbd_species;
        classes=3,
        inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
        subplot=2,
        xlab="Networks LCBD",
        ylab="Species LCBD",
        guidefontsize=7,
        bv_pal_4...
    )
    push!(biv_plots, bp)
end
titles = ["Mean", "Mean > cutoff", "Rnd", "Rnd > cutoff"]
for (bp, t) in zip(biv_plots, titles)
    plot!(bp; title=[t ""])
end
plot(biv_plots..., size = (900, 600))
savefig("lcbd_bivariate_all")

## Other maps

# Proportion of realized links
plot(links; c=:cividis, title="Proportion of realized links")
savefig(joinpath("figures", "links_proportion.png"))

# Map & compare LCBD values
plot(
    plot(lcbd_species, leg=false, c=:viridis, title="Species LCBD"),
    plot(lcbd_networks, leg=false, c=:viridis, title="Networks LCBD"),
    layout=(2,1),
    size=(600,600)
)
savefig(joinpath("figures", "lcbd_two-panels.png"))

# Univariate rescaled LCBD
plot(
    plot(rescale(lcbd_species, collect(0.0:0.05:1.0)); c=cgrad([p0, bv_pal_4[3]])),
    plot(rescale(lcbd_networks, collect(0.0:0.05:1.0)); c=cgrad([p0, bv_pal_4[2]])),
    title=["Species LCBD (rescaled)" "Networks LCBD (rescaled)"],
    layout=(2,1),
    size=(600,600)
)
savefig(joinpath("figures", "lcbd_two-panels_rescaled.png"))

## Relationship

# Visualize LCBD relationship
histogram2d(
    rescale(lcbd_networks, collect(0.0:0.05:1.0)),
    rescale(lcbd_species, collect(0.0:0.05:1.0));
    bins=20,
    xaxis=((0, 1), "Networks LCBD"),
    yaxis=((0, 1), "Species LCBD")
)
savefig(joinpath("figures", "relationship_lcbd.png"))

# Links relationship
histogram2d(
    rescale(links, collect(0.0:0.05:1.0)),
    rescale(lcbd_networks, collect(0.0:0.05:1.0));
    bins=20,
    xaxis=((0, 1), "Proportion of realized links"),
    yaxis=((0, 1), "Networks LCBD")
)
savefig(joinpath("figures", "relationship_links.png"))

# Richness relationship
histogram2d(
    rescale(Sμ, collect(0.0:0.05:1.0)),
    rescale(lcbd_networks, collect(0.0:0.05:1.0));
    bins=20,
    xaxis=((0, 1), "Proportion of realized links"),
    yaxis=((0, 1), "Networks LCBD")
)
savefig(joinpath("figures", "relationship_richness.png"))
