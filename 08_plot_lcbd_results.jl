#### Plot richness & LCBD results ####

CAN = true
include("A0_required.jl");

# Set corresponding resolution
if (@isdefined CAN) && CAN == true
    res = 2.5
else
    res = 10.0
end

# Load CairoMakie if exporting figures
if (@isdefined SAVE) && SAVE == true
    CairoMakie.activate!()
end

# Load LCBD results
include("x_load_lcbd_results.jl");

# Check results
S_all
lcbd_species_all
lcbd_networks_all

## Richness plots

# Richness for mean only
begin
    fig = background_map()
    hm2 = surface!(S_all["mean"]; colormap=:cividis, shading=false)
    Colorbar(fig[1,end+1], hm2; height=Relative(0.5), label="Expected Richness")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "richness_mean.png"), fig; px_per_unit=3.0)
end

# Richness variance for mean only
begin
    fig = background_map()
    hm2 = surface!(Sv; colormap=:cividis, shading=false)
    Colorbar(fig[1,end+1], hm2; height=Relative(0.5), label="Richness variance")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "richness_var.png"), fig; px_per_unit=3.0)
end

# Bivariate richness map
begin
    bivariate(S_all["mean"], Sv, ws; quantiles=true, classes=3, bv_pal_2...)
    bivariatelegend!(
        S_all["mean"],
        Sv;
        classes=3,
        # inset=(1, bbox(0.0, 0.17, 0.10, 0.28, :top, :right)),
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Expected richness",
        ylab="Richness variance",
        guidefontsize=7,
        bv_pal_2...
    )
end
if Makie.current_backend() == CairoMakie
    savefig(joinpath("figures", "richness_bivariate.png"))
end

## LCBD plots

# Species LCBD
begin
    fig = background_map()
    hm2 = surface!(lcbd_species_all["mean"]; colormap=:viridis, shading=false)
    Colorbar(fig[1,end+1], hm2; height=Relative(0.5), label="Relative species LCBD")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "lcbd_mean_species.png"), fig; px_per_unit=3.0)
end

# Network LCBD
begin
    fig = background_map()
    hm2 = surface!(lcbd_networks_all["mean"]; colormap=:viridis, shading=false)
    Colorbar(fig[1,end+1], hm2; height=Relative(0.5), label="Relative network LCBD")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "lcbd_mean_networks.png"), fig; px_per_unit=3.0)
end

# Bivariate species-networks LCBD for mean only
begin
    bivariate(
        lcbd_networks_all["mean"], lcbd_species_all["mean"], ws;
        quantiles=true, bv_pal_4..., classes=3,
    )
    bivariatelegend!(
        lcbd_networks_all["mean"],
        lcbd_species_all["mean"];
        classes=3,
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Networks LCBD",
        ylab="Species LCBD",
        guidefontsize=7,
        bv_pal_4...
    )
end
if Makie.current_backend() == CairoMakie
    savefig(joinpath("figures", "lcbd_bivariate_mean.png"))
end

## Sampling options

options = reshape(options, (2,2))
titles = reshape(titles, (2,2))

# All richness options
begin
    layers_all = S_all
    cbmin = mapreduce(minimum, min, values(layers_all))
    cbmax = mapreduce(maximum, max, values(layers_all))
    fig = Figure()
    for i in 1:2, j in 1:2
        o = options[i,j]
        l = layers_all[o]
        t = titles[i,j]
        hm = heatmap(
            fig[i,j], l; colormap=:cividis, colorrange=(cbmin, cbmax), axis=(;title=t)
        )
    end
    Colorbar(
        fig[:,end+1];
        height=Relative(0.5),
        colormap=:cividis,
        colorrange=(cbmin, cbmax),
        label="Expected Richness",
    )
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "sampling_options", "richness_all.png"), fig; px_per_unit=3.0)
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
    save(joinpath("figures", "sampling_options", "lcbd_species_all.png"), fig; px_per_unit=3.0)
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
    save(joinpath("figures", "sampling_options", "lcbd_networks_all.png"), fig; px_per_unit=3.0)
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
