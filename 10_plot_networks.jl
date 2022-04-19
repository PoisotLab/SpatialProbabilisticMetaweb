#### Plot network measures ####

include("A0_required.jl")

# Objects
@load joinpath("data", "jld2", "network_layers.jld2") layer layer_thr layer_rnd layer_rnd_thr
layers_all = [layer, layer_thr, layer_rnd, layer_rnd_thr]
Co = broadcast(connectance, layer)
L = broadcast(links, layer)
Lv = broadcast(links_var, layer)
Ld = broadcast(linkage_density, layer)
S = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_mean.tif"))
Sσ = geotiff(SimpleSDMPredictor, joinpath("data", "results", "richness_uncertainty.tif"))

## Some plots

# Links
plot(L; c=:acton, title="Expected number of links")
savefig(joinpath("figures", "links_mean.png"))

# Links with other color palette
plot_options = (
    cbtitle="Expected number of links", xaxis="Longitude", yaxis="Latitude", size=(650,400)
)
plot(L; c=:acton, plot_options...)
savefig(joinpath("figures", "links_mean_acton.png"))

# Link variance
plot(Lv; c=:acton, cb_title="Link variance")
savefig(joinpath("figures", "links_var.png"))

# Link bivariate map
bivariate(
    L, Lv;
    quantiles=true, classes=3, xlab="Longitude", ylab="Latitude", bv_pal_2...
)
bivariatelegend!(
    L,
    Lv;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Links",
    ylab="Link variance",
    guidefontsize=7,
    bv_pal_2...
)
plot!(title=["Links & uncertainty bivariate" ""])
savefig(joinpath("figures", "links_bivariate.png"))

# Links relationship
histogram2d(
    L,
    Lv;
    bins=20,
    xaxis=("Links"),
    yaxis=("Link variance")
)
savefig(joinpath("figures", "links_relationship.png"))

# Link coefficient of variation
Lcv = sqrt(Lv)/L
plot(Lcv; c=:cividis, title="Link coefficient of variation")
savefig(joinpath("figures", "links_coeff_var.png"))

# Link inverse-coefficient of variation or signal-to-noise ratio (SNR)
Lsnr = L/sqrt(Lv)
plot(Lsnr; c=:cividis, title="Link signal-to-noise ratio")
savefig(joinpath("figures", "links_coeff_var_inv.png"))

# Links coefficient of variation relationship
histogram2d(
    L,
    Lcv;
    bins=20,
    xaxis=("Links"),
    yaxis=("Link coefficient of variation")
)

# Link bivariate map
bivariate(
    L, Lcv;
    quantiles=true, classes=3, xlab="Longitude", ylab="Latitude", bv_pal_2...
)
bivariatelegend!(
    L,
    Lcv;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Links",
    ylab="Link coefficient of variation",
    guidefontsize=7,
    bv_pal_2...
)

## Richness

# Richness-link relationship
histogram2d(S, L, xlab="Richness", ylab="Links")
scatter(S, L, xlab="Richness", ylab="Links", alpha=0.1, legend=:none, c=:black)
savefig(joinpath("figures", "richness_relationship.png"))
plot!(xaxis=("Richness (log)", :log), yaxis=("Links (log)", :log), c=:black)
savefig(joinpath("figures", "richness_relationship_log.png"))

# Richness-link bivariate map
bivariate(
    S, L;
    quantiles=true, classes=3, xlab="Longitude", ylab="Latitude", bv_pal_2...
)
bivariatelegend!(
    S,
    L;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Richness",
    ylab="Links",
    guidefontsize=7,
    bv_pal_2...
)
savefig(joinpath("figures", "bivariate_richness_links.png"))

# Richness-link uncertainty bivariate map
bivariate(
    broadcast(v -> v^2, Sσ), Lv;
    quantiles=true, classes=3, xlab="Longitude", ylab="Latitude", bv_pal_2...
)
bivariatelegend!(
    broadcast(v -> v^2, Sσ),
    Lv;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Richness variance",
    ylab="Link variance",
    guidefontsize=7,
    bv_pal_2...
)
savefig(joinpath("figures", "bivariate_richness_links_variance.png"))

# Richness coefficient of variation
Scv = Sσ/S
plot(Scv; c=:cividis, title="Richness coefficient of variation")

# Richness-link coefficient of variation bivariate map
bivariate(
    Scv, Lcv;
    quantiles=true, classes=3, xlab="Longitude", ylab="Latitude", bv_pal_2...
)
bivariatelegend!(
    Scv,
    Lcv;
    classes=3,
    inset=(1, bbox(0.04, 0.05, 0.28, 0.28, :top, :right)),
    subplot=2,
    xlab="Richness coefficient of variation",
    ylab="Link coefficient of variation",
    guidefontsize=6,
    bv_pal_2...
)
savefig(joinpath("figures", "bivariate_richness_links_coeff.png"))

## LCBD & network measures

# Load LCBD results
include("x_load_results.jl")

# LCBD-richness relationships
begin
    scatter(S, lcbd_species_all["mean"], alpha=0.2, label="Species LCBD")
    scatter!(S, lcbd_networks_all["mean"], alpha=0.2, label="Network LCBD")
    plot!(xaxis=("Richness (log)", :log), yaxis=("LCBD"), legend=:bottomright)
end
savefig(joinpath("figures", "lcbd_relationship_richness.png"))

# LCBD-links relationships
begin
    scatter(L, lcbd_species_all["mean"], alpha=0.2, label="Species LCBD")
    scatter!(L, lcbd_networks_all["mean"], alpha=0.2, label="Network LCBD")
    plot!(xaxis=("Links (log)", :log), yaxis=("LCBD"), legend=:bottomright)
end
savefig(joinpath("figures", "lcbd_relationship_links.png"))

# LCBD-connectance relationships
begin
    scatter(Co, lcbd_species_all["mean"], alpha=0.2, label="Species LCBD")
    scatter!(Co, lcbd_networks_all["mean"], alpha=0.2, label="Network LCBD")
    plot!(xaxis=("Connectance (log)", :log), yaxis=("LCBD"), legend=:bottomright)
end
savefig(joinpath("figures", "lcbd_relationship_connectance.png"))

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
    plot(xaxis=("Richness (log)", :log), yaxis=("Links (log)", :log), legend=:bottomright)
    scatter!(S[keys(sites_mid)], L[keys(sites_mid)], label="Middle sites", alpha=0.1, c=:black)
    scatter!(S[keys(sites3)], L[keys(sites3)], label="Unique species",alpha=0.2, c=biv_colors[3])
    scatter!(S[keys(sites7)], L[keys(sites7)], label="Unique networks",alpha=0.2, c=biv_colors[7])
end
savefig(joinpath("figures", "lcbd_bivariate_scatter.png"))

# Density comparison for richness
plot(xlab="Richness", ylab="Density")
density!(S[keys(sites_mid)], label="Middle sites", c=:black)
density!(S[keys(sites3)], label="Unique species", c=biv_colors[3])
density!(S[keys(sites7)], label="Unique networks", c=biv_colors[7])
savefig(joinpath("figures", "lcbd_bivariate_density_richness.png"))

# Density comparison for links
plot(xlab="Links", ylab="Density")
density!(L[keys(sites_mid)], label="Middle sites", c=:black)
density!(L[keys(sites3)], label="Unique species", c=biv_colors[3])
density!(L[keys(sites7)], label="Unique networks", c=biv_colors[7])
savefig(joinpath("figures", "lcbd_bivariate_density_links.png"))

# Comparison for all unique species regardless of networks
begin
    _unique_spe = union(keys(sites[3]), keys(sites[6]), keys(sites[9]))
    _non_unique_spe = setdiff(keys(S), _unique_spe)
    _p1 = plot(xaxis=("Richness (log)", :log), yaxis=("Links (log)", :log), legend=:bottomright)
    scatter!(S[_non_unique_spe], L[_non_unique_spe], label="Non unique sites", alpha=0.1, c=:black)
    scatter!(S[_unique_spe], L[_unique_spe], label="Unique species", alpha=0.1, c=biv_colors[3])
    _p2 = plot(xlab="Richness", ylab="Density")
    density!(S[_non_unique_spe], label="Non unique sites", c=:black)
    density!(S[_unique_spe], label="Unique species", c=biv_colors[3])
    plot(_p1, _p2, size=(800, 400), left_margin=3mm, bottom_margin=3mm)
end
savefig(joinpath("figures", "lcbd_bivariate_unique_species.png"))

# Comparison for all unique networks regardless of species
begin
    _unique_net = union(keys(sites[7]), keys(sites[8]), keys(sites[9]))
    _non_unique_net = setdiff(keys(S), _unique_net)
    _p1 = plot(xaxis=("Richness (log)", :log), yaxis=("Links (log)", :log), legend=:bottomright)
    scatter!(S[_non_unique_net], L[_non_unique_net], label="Non unique sites", alpha=0.1, c=:black)
    scatter!(S[_unique_net], L[_unique_net], label="Unique networks", alpha=0.1, c=biv_colors[7])
    _p2 = plot(xlab="Richness", ylab="Density")
    density!(S[_non_unique_net], label="Non unique sites", c=:black)
    density!(S[_unique_net], label="Unique networks", c=biv_colors[7])
    plot(_p1, _p2, size=(800, 400), left_margin=3mm, bottom_margin=3mm)
end
savefig(joinpath("figures", "lcbd_bivariate_unique_networks.png"))

## Compare sampling options

# Links
L_all = [broadcast(links, l) for l in layers_all]
clim1 = mapreduce(minimum, min, values(L_all))
clim2 = mapreduce(maximum, max, values(L_all))
lims = (clim1, clim2)
titles = ["Mean" "Mean > cutoff" "Rnd" "Rnd > cutoff"] # for plots later on
plot(
    [plot(L; c=:acton, clim=lims) for L in L_all]...;
    # [plot(broadcast(links, l); c=:cividis) for l in layers_all]...;
    title = titles,
    cbtitle="Links",
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "links_mean_all.png"))

# Link variance
Lv_all = [broadcast(links_var, l) for l in layers_all]
clim1 = mapreduce(minimum, min, values(Lv_all))
clim2 = mapreduce(maximum, max, values(Lv_all))
lims = (clim1, clim2)
titles = ["Mean" "Mean > cutoff" "Rnd" "Rnd > cutoff"] # for plots later on
plot(
    [plot(Lv; c=:acton, clim=lims) for Lv in Lv_all]...;
    title = titles,
    cbtitle="Link variance",
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "links_var_all.png"))

## Quebec background
spatialrange = (left=-80.0, right=-50.0, bottom=45.0, top=65.0)
plot(;
    frame=:box,
    xlim=(spatialrange.left, spatialrange.right),
    ylim=(spatialrange.bottom, spatialrange.top),
    dpi=600,
)
plot!(worldshape(50), c=:lightgrey, lc=:lightgrey, grid=:none, frame=:none)
savefig("qcbackground.png")