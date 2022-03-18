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
plot(L; c=:cividis, title="Expected number of links")
savefig(joinpath("figures", "links_mean.png"))

# Links with other color palette
plot(L; c=:cividis, xaxis="Longitude", yaxis="Latitude")
savefig(joinpath("figures", "links_mean_cividis.png"))
plot(L; c=:acton, xaxis="Longitude", yaxis="Latitude")
savefig(joinpath("figures", "links_mean_acton.png"))
plot(L; c=:tokyo, xaxis="Longitude", yaxis="Latitude")
savefig(joinpath("figures", "links_mean_tokyo.png"))
plot(L; c=:viridis, xaxis="Longitude", yaxis="Latitude")
savefig(joinpath("figures", "links_mean_viridis.png"))

# Link variance
plot(Lv; c=:cividis, title="Link variance")
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
scatter(S, L, xlab="Richness", ylab="Links", alpha=0.1, legend=:none)
savefig(joinpath("figures", "richness_relationship.png"))
plot!(xaxis=("Richness (log)", :log), yaxis=("Links (log)", :log))
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
savefig(joinpath("figures", "bivariate_richness_links_uncertainty.png"))

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

## Compare sampling options

# Links
L_all = [broadcast(links, l) for l in layers_all]
clim1 = mapreduce(minimum, min, values(L_all))
clim2 = mapreduce(maximum, max, values(L_all))
lims = (clim1, clim2)
titles = ["Mean" "Mean > cutoff" "Rnd" "Rnd > cutoff"] # for plots later on
plot(
    [plot(L; c=:cividis, clim=lims) for L in L_all]...;
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
    [plot(Lv; c=:cividis, clim=lims) for Lv in Lv_all]...;
    title = titles,
    cbtitle="Link variance",
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "links_var_all.png"))
