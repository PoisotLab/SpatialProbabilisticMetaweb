#### Plot network measures ####

CAN = true
include("A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

# Objects
Co = geotiff(SimpleSDMPredictor, joinpath(results_path, "connectance.tif"))
L = geotiff(SimpleSDMPredictor, joinpath(results_path, "links_mean.tif"))
Lv = geotiff(SimpleSDMPredictor, joinpath(results_path, "links_var.tif"))
Ld = geotiff(SimpleSDMPredictor, joinpath(results_path, "links_density.tif"))
S = geotiff(SimpleSDMPredictor, joinpath(results_path, "richness_mean.tif"))
Sσ = geotiff(SimpleSDMPredictor, joinpath(results_path, "richness_uncertainty.tif"))

# Load worldshape shapefile to use as background on maps
ws = worldshape(50)

## Some plots

# Connectance
plot(Co, ws; c=:acton, cbtitle="Connectance", size=(650,400))
savefig(joinpath("figures", "links_connectance.png"))

# Links
plot(L, ws; c=:acton, cbtitle="Expected number of links")
savefig(joinpath("figures", "links_mean.png"))

# Link variance
plot(Lv, ws; c=:acton, cb_title="Link variance")
savefig(joinpath("figures", "links_var.png"))

# Link bivariate map
begin
    bivariate(L, Lv, ws; quantiles=true, classes=3, bv_pal_2...)
    bivariatelegend!(
        L,
        Lv;
        classes=3,
        # inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Links",
        ylab="Link variance",
        guidefontsize=7,
        bv_pal_2...
    )
end
savefig(joinpath("figures", "links_bivariate.png"))

# Link coefficient of variation
Lcv = sqrt(Lv)/L
plot(Lcv, ws; c=:cividis, cbtitle="Link coefficient of variation")
savefig(joinpath("figures", "links_coeff_var.png"))

# Link inverse-coefficient of variation or signal-to-noise ratio (SNR)
Lsnr = L/sqrt(Lv)
plot(Lsnr, ws; c=:cividis, cbtitle="Link signal-to-noise ratio")
savefig(joinpath("figures", "links_coeff_var_inv.png"))

## Richness

# Richness-link relationship
histogram2d(S, L, xlab="Richness", ylab="Links", cbtitle="Sites")
savefig(joinpath("figures", "richness_relationship.png"))
histogram2d(S, L; xaxis=("Richness (log)", :log), yaxis=("Links (log)", :log), cbtitle="Sites")
savefig(joinpath("figures", "richness_relationship_log.png"))

# Richness-link bivariate map
begin
    bivariate(S, L, ws; quantiles=true, classes=3, bv_pal_2...)
    bivariatelegend!(
        S,
        L;
        classes=3,
        # inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Richness",
        ylab="Links",
        guidefontsize=7,
        bv_pal_2...
    )
end
savefig(joinpath("figures", "bivariate_richness_links.png"))

# Richness-link uncertainty bivariate map
begin
    bivariate(broadcast(v -> v^2, Sσ), Lv, ws; quantiles=true, classes=3, bv_pal_2...)
    bivariatelegend!(
        broadcast(v -> v^2, Sσ),
        Lv;
        classes=3,
        # inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Richness variance",
        ylab="Link variance",
        guidefontsize=7,
        bv_pal_2...
    )
end
savefig(joinpath("figures", "bivariate_richness_links_variance.png"))

## LCBD & network measures

# Load LCBD results
include("x_load_lcbd_results.jl");

# LCBD-richness relationships
plot(
    histogram2d(S, lcbd_species_all["mean"]; cb=:none),
    histogram2d(S, lcbd_networks_all["mean"]; cb=:none);
    xaxis="Richness", yaxis="Relative LCBD", title=["Species LCBD" "Networks LCBD"],
    size=(700, 400)
)
savefig(joinpath("figures", "lcbd_relationship_richness.png"))

# LCBD-links relationships
plot(
    histogram2d(L, lcbd_species_all["mean"]; cb=:none),
    histogram2d(L, lcbd_networks_all["mean"]; cb=:none);
    xaxis=("Links (log)", :log), yaxis="Relative LCBD", title=["Species LCBD" "Networks LCBD"],
    size=(700, 400)
)
savefig(joinpath("figures", "lcbd_relationship_links.png"))

## Compare link & richness density of unique sites

# Extract the bivariate values
biv_layer, biv_colors = get_bivariate_values(
    lcbd_networks_all["mean"],
    lcbd_species_all["mean"];
    bv_pal_4...
)
# plot(convert(Float32, biv_layer), c=biv_colors)

# Get the specific sites for each group
sites3 = broadcast(v -> v == 3 ? 1 : nothing, biv_layer)
sites7 = broadcast(v -> v == 7 ? 1 : nothing, biv_layer)
sites_mid = broadcast(v -> v != 3 && v != 7 ? 1 : nothing, biv_layer)
sites = [broadcast(v -> v == i ? true : nothing, biv_layer) for i in 1:9]
union(keys(sites[3]), keys(sites[6]), keys(sites[9]))

# Plot the two extremas in a different color
begin
    _p1 = plot(xaxis=("Richness (log)", :log), yaxis=("Links (log)", :log), legend=:bottomright)
    scatter!(S[keys(S)], L[keys(S)], label="Other sites", alpha=0.1, c=:black)
    scatter!(S[keys(sites3)], L[keys(sites3)], label="Unique species only",alpha=0.2, c=biv_colors[3])
    scatter!(S[keys(sites7)], L[keys(sites7)], label="Unique networks only",alpha=0.2, c=biv_colors[7])
    _p2 = plot(xaxis=("Richness"), yaxis=("Density"))
    density!(S[keys(sites_mid)], label="Other sites", c=:black)
    density!(S[keys(sites3)], label="Unique species only", c=biv_colors[3])
    density!(S[keys(sites7)], label="Unique networks only", c=biv_colors[7])
    _p3 = plot(xaxis=("Links"), yaxis=("Density"))
    density!(L[keys(sites_mid)], label="Other sites", c=:black)
    density!(L[keys(sites3)], label="Unique species only", c=biv_colors[3])
    density!(L[keys(sites7)], label="Unique networks only", c=biv_colors[7])
    _layout = @layout [a [b; c]]
    plot(_p1, _p2, _p3; size=(800, 400), left_margin=3mm, bottom_margin=3mm, layout=_layout)
end; # do not display as VS Code might crash
savefig(joinpath("figures", "lcbd_bivariate_densities.png"))

_layout = @layout [a [b; c]]
plot(
    plot(), plot(), plot(); layout=_layout
)

## Compare sampling options

# layers_all = [layer, layer_thr, layer_rnd, layer_rnd_thr]
# L_all = [broadcast(links, l) for l in layers_all]

# Links for all options
# clim1 = mapreduce(minimum, min, values(L_all))
# clim2 = mapreduce(maximum, max, values(L_all))
# L_all_plots = []
# lims = (clim1, clim2)
# plot(
#     [plot(L, ws; c=:acton, clim=lims) for L in L_all]...;
#     cbtitle="Expected number of links",
#     layout=(2,2),
#     size=(1000,600),
#     left_margin=3mm
# )
# savefig(joinpath("figures", "sampling_options", "links_all.png"))

# Link variance
# Lv_all = [broadcast(links_var, l) for l in layers_all]
# clim1 = mapreduce(minimum, min, values(Lv_all))
# clim2 = mapreduce(maximum, max, values(Lv_all))
# lims = (clim1, clim2)
# titles = ["Mean" "Mean > cutoff" "Rnd" "Rnd > cutoff"] # for plots later on
# plot(
#     [plot(Lv, ws; c=:acton, clim=lims) for Lv in Lv_all]...;
#     title = titles,
#     cbtitle="Link variance",
#     layout=(2,2),
#     size=(900,600),
# )
# savefig(joinpath("figures", "sampling_options", "links_var_all.png"))
