#### Plot richness & LCBD results ####

include("A0_required.jl")

## Load data

# Define names for the 4 sampling options
options = ["mean", "mean_thr", "rand", "rand_thr"]
titles = ["Mean" "Mean > cutoff" "Rnd" "Rnd > cutoff"] # for plots later on

# Richness layers
S_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath("data", "results", "richness_$(opt).tif")
    S_all[opt] = geotiff(SimpleSDMPredictor, path)
end
S_all

# Species LCBD layers
lcbd_species_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath("data", "results", "lcbd_species_$(opt).tif")
    lcbd_species_all[opt] = geotiff(SimpleSDMPredictor, path)
end
lcbd_species_all

# Networks LCBD layers
lcbd_networks_all = Dict{String, SimpleSDMPredictor}()
for opt in options
    path = joinpath("data", "results", "lcbd_networks_$(opt).tif")
    lcbd_networks_all[opt] = geotiff(SimpleSDMPredictor, path)
end
lcbd_networks_all

# Others
Sσ = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_uncertainty.tif"))
spatialrange = (left=-80.0, right=-50.0, bottom=45.0, top=65.)
reference_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; spatialrange...)

## Richness plots

# All richness options
clim1 = mapreduce(minimum, min, values(S_all))
clim2 = mapreduce(maximum, max, values(S_all))
lims = (clim1, clim2)
plot(
    [plot(S_all[opt]; c=:cividis, clim=lims) for opt in options]...;
    title = titles,
    cbtitle="Species richness",
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "richness_all.png"))
# For committee document
plot(
    [plot(S_all[opt]; c=:cividis, clim=lims) for opt in options]...;
    cbtitle="Species richness",
    layout=(2,2),
    size=(1000,600),
    xaxis="Latitude",
    yaxis="Longitude",
    left_margin=3mm
)
savefig(joinpath("figures", "richness_all_committee.png"))
plot!(title = ["a)" "b)" "c)" "d)"], titlepos=:left)
savefig(joinpath("figures", "richness_all_committee_lettered.png"))

# Richness for mean only
plot(S_all["mean"], c=:cividis, cbtitle="Expected richness", size=(650, 400))
plot!(xaxis="Longitude", yaxis="Latitude")
savefig(joinpath("figures", "richness_mean.png"))

# Richness variance for mean only
plot(Sσ, c=:cividis, cbtitle="Richness variance", size=(650, 400))
plot!(xaxis="Longitude", yaxis="Latitude")
savefig(joinpath("figures", "richness_var.png"))

# Univariate richness maps
plot(
    plot(S_all["mean"], title="Expected richness", c=cgrad([p0, bv_pal_2[2]])),
    plot(Sσ, title="Std. dev. of richness", c=cgrad([p0, bv_pal_2[3]]));
    layout=(2,1),
    size=(600, 600)
)
savefig(joinpath("figures", "richness_two-panels.png"))

# Bivariate richness map
bivariate(
    S_all["mean"], Sσ;
    quantiles=true, classes=3, xlab="Longitude", ylab="Latitude", bv_pal_2...
)
bivariatelegend!(
    S_all["mean"],
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
    [plot(lcbd_species_all[opt]/maximum(lcbd_species_all[opt]); c=:viridis) for opt in options]...;
    title=titles,
    cbtitle="Relative species LCBD",
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "lcbd_species_all.png"))

# Networks LCBD plots

# All networks LCBD options
netw_plots = []
for opt in options
    _l = lcbd_networks_all[opt]
    _lrel = _l/maximum(_l)
    p = plot(reference_layer; c=:lightgrey)
    plot!(p, _lrel; c=:viridis, clim=extrema(_lrel))
    push!(netw_plots, p)
end
plot(netw_plots..., size = (900, 600), title=titles, cbtitle="Relative networks LCBD")
savefig(joinpath("figures", "lcbd_networks_all.png"))

# Bivariate species-networks LCBD
biv_plots = []
for (i, opt) in enumerate(options)
    bp = bivariate(
        lcbd_networks_all[opt], lcbd_species_all[opt];
        quantiles=true, bv_pal_4..., classes=3, title=titles[i]
    )
    bp = bivariatelegend!(
        lcbd_networks_all[opt],
        lcbd_species_all[opt];
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
plot(biv_plots..., size = (900, 600))
savefig(joinpath("figures", "lcbd_bivariate_all.png"))

# Bivariate species-networks LCBD for mean only
bivariate(
    lcbd_networks_all["mean"], lcbd_species_all["mean"];
    quantiles=true, bv_pal_4..., classes=3,
    xaxis="Longitude", yaxis="Latitude"
)
bivariatelegend!(
    lcbd_networks_all["mean"],
    lcbd_species_all["mean"];
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Networks LCBD",
    ylab="Species LCBD",
    guidefontsize=7,
    bv_pal_4...
)
savefig(joinpath("figures", "lcbd_bivariate_mean.png"))

## Other maps

# Map & compare LCBD values
plot(
    plot(lcbd_species_all["mean"], leg=false, c=:viridis, title="Species LCBD"),
    plot(lcbd_networks_all["mean"], leg=false, c=:viridis, title="Networks LCBD"),
    layout=(2,1),
    size=(600,600)
)
savefig(joinpath("figures", "lcbd_two-panels.png"))

# Plot separately
plot_options = (cb=:none, xticks=:none, yticks=:none, frame=:box)
plot(lcbd_species_all["mean"]; c=:viridis, plot_options...)
savefig(joinpath("figures", "lcbd_mean_species.png"))
plot(lcbd_networks_all["mean"]; c=:viridis, plot_options...)
savefig(joinpath("figures", "lcbd_mean_networks.png"))

# Univariate rescaled LCBD
plot(
    plot(rescale(lcbd_species_all["mean"], collect(0.0:0.05:1.0)); c=cgrad([p0, bv_pal_4[3]])),
    plot(rescale(lcbd_networks_all["mean"], collect(0.0:0.05:1.0)); c=cgrad([p0, bv_pal_4[2]])),
    title=["Species LCBD (rescaled)" "Networks LCBD (rescaled)"],
    layout=(2,1),
    size=(600,600)
)
savefig(joinpath("figures", "lcbd_two-panels_rescaled.png"))

# Plot separately
plot_options = (xaxis="Longitude", yaxis="Latitude", size=(650,400))
plot(rescale(lcbd_species_all["mean"], collect(0.0:0.05:1.0)); c=cgrad([p0, bv_pal_4[3]]))
plot!(; cbtitle= "Species LCBD (rescaled)", plot_options...)
savefig(joinpath("figures", "lcbd_mean_rescaled_species.png"))
plot(rescale(lcbd_networks_all["mean"], collect(0.0:0.05:1.0)); c=cgrad([p0, bv_pal_4[2]]))
plot!(; cbtitle= "Network LCBD (rescaled)", plot_options...)
savefig(joinpath("figures", "lcbd_mean_rescaled_networks.png"))

## Relationship

# Visualize LCBD relationship
histogram2d(
    rescale(lcbd_networks_all["mean"], collect(0.0:0.05:1.0)),
    rescale(lcbd_species_all["mean"], collect(0.0:0.05:1.0));
    bins=20,
    xaxis=((0, 1), "Networks LCBD"),
    yaxis=((0, 1), "Species LCBD")
)
savefig(joinpath("figures", "relationship_lcbd.png"))

# Richness relationship
histogram2d(
    rescale(S_all["mean"], collect(0.0:0.05:1.0)),
    rescale(lcbd_networks_all["mean"], collect(0.0:0.05:1.0));
    bins=20,
    xaxis=((0, 1), "Richness"),
    yaxis=((0, 1), "Networks LCBD")
)
savefig(joinpath("figures", "relationship_richness.png"))
