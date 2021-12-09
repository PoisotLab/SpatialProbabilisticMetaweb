#### Bivariate maps ####

include("A0_required.jl")

## Load data

# Richness layers
Sμ = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_mean.tif"))
Sσ = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_uncertainty.tif"))
Sr = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_rand.tif"))
Sμ_cut = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_mean_thr.tif"))
Sr_cut = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_rand_thr.tif"))

# Species LCBD layers
lcbd_layers = fill(SimpleSDMPredictor(rand(Float32, 2,2)), 4)
lcbd_layers[1] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_species_mean.tif"))
lcbd_layers[2] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_species_rand.tif"))
lcbd_layers[3] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_species_mean_thr.tif"))
lcbd_layers[4] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_species_rand_thr.tif"))
lcbd_species = lcbd_layers[1]

# Networks LCBD layers
lcbd_networks = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_networks_mean.tif"))
lcbd_networks_rnd = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_networks_rand.tif"))
lcbd_networks_thr = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_networks_mean_thr.tif"))
lcbd_networks_rnd_thr = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_networks_rand_thr.tif"))
L = geotiff(SimpleSDMPredictor, joinpath("data", "results", "links.tif"))

## Bivariate maps

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
savefig(joinpath("figures", "lcbd_bivariate_all.png"))

# Networks LCBD
netw_plots = []
for lcbd_n in [lcbd_networks, lcbd_networks_thr, lcbd_networks_rnd, lcbd_networks_rnd_thr]
    p = plot(lcbd_networks; c=:lightgrey)
    plot!(p, lcbd_n; c=:viridis, clim=extrema(lcbd_n))
    push!(netw_plots, p)
end
titles = ["Mean", "Mean > cutoff", "Rnd", "Rnd > cutoff"]
plot(netw_plots..., size = (900, 600), title=permutedims(titles))
savefig(joinpath("figures", "lcbd_networks_all.png"))

## Other maps

# Proportion of realized links
plot(L; c=:cividis, title="Proportion of realized links")
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
    rescale(L, collect(0.0:0.05:1.0)),
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
