## Probabilistic distributions

include("04_aggregate_sdms.jl")

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
# networks = assemble_networks(reference_layer, P, D, A, cutoffs; type="avg_thr"); # 30 sec.
# networks = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd"); # 2 min
# networks = assemble_networks(reference_layer, P, D, A, cutoffs; type="rnd_thr"); # 30 sec.

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
lcbd_networks[keys(lcbd_networks)] = LCBD(hellinger(Zfull))[1]

# Map & compare LCBD values
plot(
    plot(lcbd_species, leg=false, c=:viridis, title="Species LCBD"),
    plot(lcbd_networks, leg=false, c=:viridis, title="Networks LCBD"),
    layout=(2,1),
    size=(600,600)
)

# Prepare bivariate colors
p0 = colorant"#e8e8e8"
bv_pal_1 = (p0=p0, p1=colorant"#64acbe", p2=colorant"#c85a5a")
bv_pal_2 = (p0=p0, p1=colorant"#73ae80", p2=colorant"#6c83b5")
bv_pal_3 = (p0=p0, p1=colorant"#9972af", p2=colorant"#c8b35a")
bv_pal_4 = (p0=p0, p1=colorant"#be64ac", p2=colorant"#5ac8c8")
# Bivariate LCBD
bivariate(lcbd_networks, lcbd_species; quantiles=true, bv_pal_4..., classes=3)
bivariatelegend!(
    lcbd_networks,
    lcbd_species;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Networks LCBD",
    ylab="Species LCBD",
    guidefontsize=7,
    bv_pal_4...
)

# Univariate rescaled LCBD
plot(rescale(lcbd_networks, collect(0.0:0.05:1.0)); c=cgrad([p0, bv_pal_4[2]]))
plot(rescale(lcbd_species, collect(0.0:0.05:1.0)); c=cgrad([p0, bv_pal_4[3]]))

# Visualize relationship
histogram2d(
    rescale(lcbd_networks, collect(0.0:0.05:1.0)),
    rescale(lcbd_species, collect(0.0:0.05:1.0));
    bins=20
)
xaxis!((0,1))
yaxis!((0,1))