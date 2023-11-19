#### Plot richness & LCBD results ####

CAN = true
include("A0_required.jl");

# Load CairoMakie if exporting figures
if (@isdefined SAVE) && SAVE == true
    CairoMakie.activate!()
end

# Load LCBD results
include("scripts/x_load_lcbd_results.jl");

# Check results
S_all
lcbd_species_all
lcbd_networks_all

## Richness plots

# Richness for mean only
begin
    fig = background_map()
    sf = surface!(S_all["mean"]; colormap=:cividis, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Expected Richness")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "richness_mean.png"), fig)
end

# Richness variance for mean only
begin
    fig = background_map()
    sf = surface!(Sv; colormap=:cividis, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Richness variance")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "richness_var.png"), fig)
end

# Bivariate richness map
begin
    fig = Figure()
    g1 = fig[1:16, 1:4] = GridLayout()
    g2 = fig[2:5, end] = GridLayout()

    p1 = background_map(g1[1,1])
    sf = bivariatesurface!(p1, S_all["mean"], Sv; bv_pal_2...)

    p2 = Axis(g2[1,1]; aspect=1, xlabel = "Expected richness", ylabel = "Richness variance")
    l2 = bivariatelegend!(p2, S_all["mean"], Sv; bv_pal_2...)
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "richness_bivariate.png"), fig)
end

## LCBD plots

# Species LCBD
begin
    fig = background_map()
    sf = surface!(lcbd_species_all["mean"]; colormap=:viridis, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Relative species LCBD")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "lcbd_mean_species.png"), fig)
end

# Network LCBD
begin
    fig = background_map()
    sf = surface!(lcbd_networks_all["mean"]; colormap=:viridis, shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Relative network LCBD")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "lcbd_mean_networks.png"), fig)
end

# Bivariate species-networks LCBD for mean only
begin
    fig = Figure()
    g1 = fig[1:16, 1:4] = GridLayout()
    g2 = fig[2:5, end] = GridLayout()

    p1 = background_map(g1[1,1])
    sf = bivariatesurface!(p1, lcbd_species_nan, lcbd_networks_all["mean"])

    p2 = Axis(g2[1,1]; aspect = 1, xlabel = "Species LCBD", ylabel = "Network LCBD")
    l2 = bivariatelegend!(p2, lcbd_species_nan, lcbd_networks_all["mean"])
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "lcbd_bivariate_mean.png"), fig)
end

## Sampling options

options = reshape(options, (2,2))
# titles = reshape(titles, (2,2))
titles = [
    "A) Mean value" "C) Random value";
    "B) Mean value + threshold" "D) Random value + threshold"
]

# All richness options
begin
    layers_all = S_all
    cbmin = mapreduce(minimum, min, values(layers_all))
    cbmax = mapreduce(maximum, max, values(layers_all))
    fig = Figure(; resolution=(1250,600))
    for i in 1:2, j in 1:2
        o = options[i,j]
        l = layers_all[o]
        t = titles[i,j]
        p = background_map(fig[i,j]; title=t, titlealign=:left)
        s = surface!(
            fig[i,j], l; colormap=:cividis, colorrange=(cbmin, cbmax), shading=false
        )
        Colorbar(p[1,2], s; height=Relative(0.5), label="Expected Richness")
    end
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "sampling_options", "richness_all.png"), fig)
end

# All species LCBD options
#= # returns NaNs
begin
    layers_all = lcbd_species_all
    cbmin = mapreduce(minimum, min, values(layers_all))
    cbmax = mapreduce(maximum, max, values(layers_all))
    fig = Figure()
    for i in 1:2, j in 1:2
        o = options[i,j]
        l = layers_all[o]
        t = titles[i,j]
        hm = heatmap(fig[i,j], l; colorrange=(cbmin, cbmax), axis=(;title=t))
    end
    Colorbar(
        fig[:,end+1];
        height=Relative(0.5),
        label="Relative species LCBD",
        colorrange=(cbmin, cbmax)
    )
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "sampling_options", "lcbd_species_all.png"), fig)
end
=#

# All networks LCBD options
# NO DATA FOR NOW
#=
begin
    layers_all = lcbd_networks_all
    cbmin = mapreduce(minimum, min, values(layers_all))
    cbmax = mapreduce(maximum, max, values(layers_all))
    fig = Figure()
    for i in 1:2, j in 1:2
        o = options[i,j]
        l = layers_all[o]
        t = titles[i,j]
        hm = heatmap(fig[i,j], l; colorrange=(cbmin, cbmax), axis=(;title=t))
    end
    Colorbar(
        fig[:,end+1];
        height=Relative(0.5),
        label="Relative species LCBD",
        colorrange=(cbmin, cbmax)
    )
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "sampling_options", "lcbd_networks_all.png"), fig)
end

# Bivariate species-networks LCBD
biv_plots = []
for (i, opt) in enumerate(options)
    bp = bivariate(
        lcbd_networks_all[opt], lcbd_species_all[opt], ws;
        quantiles=true, bv_pal_4..., classes=3, title=titles[i]
    )
    bp = bivariatelegend!(
        lcbd_networks_all[opt],
        lcbd_species_all[opt];
        classes=3,
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Networks LCBD",
        ylab="Species LCBD",
        guidefontsize=7,
        bv_pal_4...
    )
    push!(biv_plots, bp)
end
plot(biv_plots..., size = (900, 600))
if Makie.current_backend() == CairoMakie
    savefig(joinpath("figures", "sampling_options", "lcbd_bivariate_all.png"))
end
=#
