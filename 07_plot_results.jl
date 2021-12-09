#### Bivariate maps ####

include("A0_required.jl")

## Load data

# Richness layers
S_all = fill(SimpleSDMPredictor(rand(Float64, 2,2)), 4)
S_all[1] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_mean.tif"))
S_all[2] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_rand.tif"))
S_all[3] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_mean_thr.tif"))
S_all[4] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_rand_thr.tif"))
Sσ = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_uncertainty.tif"))

# Species LCBD layers
lcbd_species_all = fill(SimpleSDMPredictor(rand(Float32, 2,2)), 4)
lcbd_species_all[1] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_species_mean.tif"))
lcbd_species_all[2] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_species_mean_thr.tif"))
lcbd_species_all[3] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_species_rand.tif"))
lcbd_species_all[4] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_species_rand_thr.tif"))
lcbd_species = lcbd_layers[1]

# Networks LCBD layers
lcbd_networks_all = fill(SimpleSDMPredictor(rand(Float32, 2,2)), 4)
lcbd_networks_all[1] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_networks_mean.tif"))
lcbd_networks_all[2] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_networks_mean_thr.tif"))
lcbd_networks_all[3] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_networks_rand.tif"))
lcbd_networks_all[4] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "lcbd_networks_rand_thr.tif"))

# Others
L = geotiff(SimpleSDMPredictor, joinpath("data", "results", "links.tif"))

## Richness plots

# All richness options
subtitles = ["Mean" "Mean > cutoff" "Rnd" "Rnd > cutoff"]
clim1 = mapreduce(minimum, min, S_all)
clim2 = mapreduce(maximum, max, S_all)
lims = (clim1, clim2)
plot(
    [plot(S; c=:cividis, title=t, clim=lims) for (S,t) in zip(S_all, subtitles)]...;
    cb_title="Species richness",
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "richness_all.png"))

# Univariate richness maps
plot(
    plot(Sμ, title="Expected richness", c=cgrad([p0, bv_pal_2[2]])),
    plot(Sσ, title="Std. dev. of richness", c=cgrad([p0, bv_pal_2[3]]));
    layout=(2,1),
    size=(600, 600)
)
savefig(joinpath("figures", "richness_two-panels.png"))

# Bivariate richness map
bivariate(Sμ, Sσ; quantiles=true, classes=3, xlab="Longitude", ylab="Latitude", bv_pal_2...)
bivariatelegend!(
    Sμ,
    Sσ;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Expected richness",
    ylab="Std. dev. of richness",
    guidefontsize=7,
    bv_pal_2...
)
plot!(title=["Richness & uncertainty bivariate" ""])
savefig(joinpath("figures", "richness_bivariate.png"))

## Species LCBD plots

# All LCBD options
plot(
    plot(lcbd_species_all[1]; c=:viridis, title="LCBD means"),
    plot(lcbd_species_all[2]; c=:viridis, title="LCBD means cut"),
    plot(lcbd_species_all[3]; c=:viridis, title="LCBD rands"),
    plot(lcbd_species_all[4]; c=:viridis, title="LCBD rands cut"),
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "lcbd_species_all.png"))

# Networks LCBD plots

# All networks LCBD options
netw_plots = []
for lcbd_n in [lcbd_networks, lcbd_networks_thr, lcbd_networks_rnd, lcbd_networks_rnd_thr]
    p = plot(lcbd_networks; c=:lightgrey)
    plot!(p, lcbd_n; c=:viridis, clim=extrema(lcbd_n))
    push!(netw_plots, p)
end
titles = ["Mean", "Mean > cutoff", "Rnd", "Rnd > cutoff"]
plot(netw_plots..., size = (900, 600), title=permutedims(titles))
savefig(joinpath("figures", "lcbd_networks_all.png"))

# Bivariate species-networks LCBD
biv_plots = []
for lcbd_n in [lcbd_networks, lcbd_networks_thr, lcbd_networks_rnd, lcbd_networks_rnd_thr]
    bp = bivariate(lcbd_n, lcbd_species[1]; quantiles=true, bv_pal_4..., classes=3)
    bp = bivariatelegend!(
        lcbd_n,
        lcbd_species[1];
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

## Other maps

# Proportion of realized links
plot(L; c=:cividis, title="Proportion of realized links")
savefig(joinpath("figures", "links_proportion.png"))

# Map & compare LCBD values
plot(
    plot(lcbd_species_all[1], leg=false, c=:viridis, title="Species LCBD"),
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
