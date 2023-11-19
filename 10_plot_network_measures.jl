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
    save(joinpath("figures", "links_connectance.png"), fig)
end

# Links
begin
    fig = background_map()
    sf = surface!(L; colormap=:acton, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Expected number of links")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "links_mean.png"), fig)
end

# Link variance
begin
    fig = background_map()
    sf = surface!(Lv; colormap=:acton, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Link variance")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "links_var.png"), fig)
end

# Link bivariate map
begin
    fig = Figure()
    g1 = fig[1:16, 1:4] = GridLayout()
    g2 = fig[2:5, end] = GridLayout()

    p1 = background_map(g1[1,1])
    sf = bivariatesurface!(p1, L, Lv; bv_pal_2...)

    p2 = Axis(g2[1,1]; aspect = 1, xlabel = "Links", ylabel = "Link variance")
    l2 = bivariatelegend!(p2, L, Lv; bv_pal_2...)
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "links_bivariate.png"), fig)
end

## Richness

# Richness-link bivariate map
begin
    fig = Figure()
    g1 = fig[1:16, 1:4] = GridLayout()
    g2 = fig[2:5, end] = GridLayout()

    p1 = background_map(g1[1,1])
    sf = bivariatesurface!(p1, S, L; bv_pal_2...)

    p2 = Axis(g2[1,1]; aspect = 1, xlabel = "Richness", ylabel = "Links")
    l2 = bivariatelegend!(p2, S, L; bv_pal_2...)
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "bivariate_richness_links.png"), fig)
end

# Richness-link uncertainty bivariate map
begin
    fig = Figure()
    g1 = fig[1:16, 1:4] = GridLayout()
    g2 = fig[2:5, end] = GridLayout()

    p1 = background_map(g1[1,1])
    sf = bivariatesurface!(p1, Sv, Lv; bv_pal_2...)

    p2 = Axis(g2[1,1]; aspect = 1, xlabel = "Richness variance", ylabel = "Link variance")
    l2 = bivariatelegend!(p2, Sv, Lv; bv_pal_2...)
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "bivariate_richness_links_variance.png"), fig)
end

## LCBD & network measures

# Load LCBD results
include(joinpath("scripts", "x_load_lcbd_results.jl"));

## Compare link & richness density of unique sites

# Extract the bivariate values
biv_layer = bivariatelayer(lcbd_species_nan, lcbd_networks_all["mean"])
biv_colors = _get_bivariate_colormap()

# Get the specific sites for each group
sites3 = findall(==(3), biv_layer)
sites7 = findall(==(7), biv_layer)
sites_mid = setdiff(keys(biv_layer), union(sites3, sites7))

# Plot the two extremas in a different color
_S3 = try Float64.(S[sites3]) catch; Vector{Float64}[] end
_S7 = Float64.(S[sites7])
_Smid = Float64.(S[sites_mid])
_L3 = try Float64.(L[sites3]) catch; Vector{Float64}[] end
_L7 = Float64.(L[sites7])
_Lmid = Float64.(L[sites_mid])
begin
    fig = Figure()
    ax1 = Axis(
        fig[1:2,1]; xlabel="Richness", ylabel="Links",
        xscale=Makie.pseudolog10, yscale=Makie.pseudolog10
    )
    scatter!(_Smid, _Lmid, label="Other sites", color=(:black, 0.1))
    try scatter!(_S3, _L3, label="Unique species only", color=(biv_colors[3], 0.2)) catch; end
    scatter!(_S7, _L7, label="Unique networks only", color=(biv_colors[7], 0.2))
    axislegend(ax1, position = :rb)
    ax2 = Axis(fig[1,2]; xlabel="Richness", ylabel="Density")
    density!(ax2, _Smid; color=(:black, 0.3), strokecolor=:black, strokewidth=3)
    try density!(ax2, _S3; color=(biv_colors[3], 0.3), strokecolor=biv_colors[3], strokewidth=3) catch; end
    density!(ax2, _S7; color=(biv_colors[7], 0.3), strokecolor=biv_colors[7], strokewidth=3)
    labels = ["Other sites", "Unique species only", "Unique networks only"]
    elements = [
        PolyElement(polycolor = (col, 0.3), polystrokecolor=col, polystrokewidth=3)
        for col in [:black, biv_colors[3], biv_colors[7]]
    ]
    axislegend(ax2, elements, labels)
    ax3 = Axis(fig[2,2]; xlabel="Links", ylabel="Density")
    density!(ax3, _Lmid; color=(:black, 0.3), strokecolor=:black, strokewidth=3)
    try density!(ax3, _L3; color=(biv_colors[3], 0.3), strokecolor=biv_colors[3], strokewidth=3) catch; end
    density!(ax3, _L7; color=(biv_colors[7], 0.3), strokecolor=biv_colors[7], strokewidth=3)
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "lcbd_bivariate_densities.png"), fig)
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
