#### Plot network measures ####

CAN = true
include("A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

# Load CairoMakie if exporting figures
if (@isdefined SAVE) && SAVE == true
    CairoMakie.activate!()
end

# Objects
Co = read_geotiff(joinpath(results_path, "connectance.tif"), SimpleSDMPredictor)
L = read_geotiff(joinpath(results_path, "links_mean.tif"), SimpleSDMPredictor)
Lv = read_geotiff(joinpath(results_path, "links_var.tif"), SimpleSDMPredictor)
Ld = read_geotiff(joinpath(results_path, "links_density.tif"), SimpleSDMPredictor)
S = read_geotiff(joinpath(results_path, "richness_mean.tif"), SimpleSDMPredictor)
Sv = read_geotiff(joinpath(results_path, "richness_uncertainty.tif"), SimpleSDMPredictor)

## Some plots

# Connectance
begin
    fig = background_map()
    sf = surface!(Co; colormap=:acton, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Connectance")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "links_connectance.png"), fig; px_per_unit=3.0)
end

# Links
begin
    fig = background_map()
    sf = surface!(L; colormap=:acton, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Expected number of links")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "links_mean.png"), fig; px_per_unit=3.0)
end

# Link variance
begin
    fig = background_map()
    sf = surface!(Lv; colormap=:acton, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Link variance")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "links_var.png"), fig; px_per_unit=3.0)
end

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
if Makie.current_backend() == CairoMakie
    savefig(joinpath("figures", "links_bivariate.png"))
end

## Richness

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
if Makie.current_backend() == CairoMakie
    savefig(joinpath("figures", "bivariate_richness_links.png"))
end

# Richness-link uncertainty bivariate map
begin
    bivariate(Sv, Lv, ws; quantiles=true, classes=3, bv_pal_2...)
    bivariatelegend!(
        Sv,
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
if Makie.current_backend() == CairoMakie
    savefig(joinpath("figures", "bivariate_richness_links_variance.png"))
end

## LCBD & network measures

# Load LCBD results
include("x_load_lcbd_results.jl");

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
if Makie.current_backend() == CairoMakie
    savefig(joinpath("figures", "lcbd_bivariate_densities.png"))
end

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
