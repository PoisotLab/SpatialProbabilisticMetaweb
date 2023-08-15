#### Ecoregion plots

CAN = true
include("A0_required.jl");

# Load the corresponding results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ecoresults_path = joinpath("data", "results", "ecoregions");
else
    ecoresults_path = joinpath("xtras", "results", "ecoregions");
end

# Load CairoMakie if exporting figures
if (@isdefined SAVE) && SAVE == true
    CairoMakie.activate!()
end

## Basic network measures

fig_path = joinpath("figures", "ecoregions")
isdir(fig_path) || mkdir(fig_path)

# Define the network measures to use
network_measures = ["Co", "L", "Lv", "Ld"]
measures = [network_measures..., "S", "Sv", "LCBD_species", "LCBD_networks"]
measures_ts = [
    "Connectance", "Number of links", "Link variance", "Linkage density",
    "Richness", "Richness variance", "Relative species LCBD", "Relative network LCBD"
]
summary_fs = ["median", "iqr89"]
summary_ts = ["median", "89% IQR"]

# Predefine set of options
opt = []
for m in measures, fs in summary_fs
    o = (m = m, fs = fs)
    push!(opt, o)
end
opt

# Load the ecoregion summary layers
ecoregion_layers = Dict{String, SimpleSDMResponse}()
for o in opt
    # Load layer
    path = joinpath(ecoresults_path, "ecoregion_$(o.m)_$(o.fs).tif")
    ecoregion_layers["$(o.m)_$(o.fs)"] = read_geotiff(path, SimpleSDMResponse)
    # Replace zero values (sites not in an ecoregion)
    replace!(ecoregion_layers["$(o.m)_$(o.fs)"], 0.0 => nothing)
end
ecoregion_layers

## Make some plots!!

# Plot results
begin
    fig = Figure(; resolution=(1200,600))
    for i in 1:2, j in 1:2
        m = reshape(network_measures, (2,2))[i,j]
        t = reshape(measures_ts[1:4], (2,2))[i,j]
        p = background_map(fig[i,j]; title=t, titlealign=:left)
        s = surface!(ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        Colorbar(p[1,2], s; height=Relative(0.5))
    end
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_all_median.png"), fig; px_per_unit=3.0)
end

# Some variations
ecoregion_plots = Dict{String, Figure}()
for (m,t) in zip(measures, measures_ts)
    begin
        f = Figure(; resolution=(850,800))

        p1 = background_map(f[1,1]; title="Median", titlealign=:left)
        p2 = background_map(f[2,1]; title="89% IQR", titlealign=:left)

        sf1 = surface!(p1, ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        sf2 = surface!(p2, ecoregion_layers["$(m)_iqr89"]; colormap=:inferno, shading=false)

        Colorbar(p1[1,2], sf1; height=Relative(0.5), label="$(t)")
        Colorbar(p2[1,2], sf2; height=Relative(0.5), label="$(t) 89% IQR")

        ecoregion_plots[m] = f;
    end;
end
ecoregion_plots["Co"]
ecoregion_plots["L"]
ecoregion_plots["Lv"]
ecoregion_plots["Ld"]
ecoregion_plots["S"]
ecoregion_plots["Sv"]
ecoregion_plots["LCBD_species"]
ecoregion_plots["LCBD_networks"]

# Export
if Makie.current_backend() == CairoMakie
    @threads for m in String.(keys(ecoregion_plots))
        save(joinpath(fig_path, "ecoregion_$m.png"), ecoregion_plots[m]; px_per_unit=3.0)
    end
end

## Compare with richness

# Compare with variance measures for the ecoregion
begin
    ms = ["S" "Sv"; "L" "Lv"]
    ts = ["Richness" "Richness variance"; "Links" "Link variance"]
    fig = Figure(; resolution=(1275,600))
    for i in 1:2, j in 1:2
        m = ms[i,j]
        t = ts[i,j]
        p = background_map(fig[i,j]; title=t, titlealign=:left)
        s = surface!(ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        Colorbar(p[1,2], s; height=Relative(0.5), label=t)
    end
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_comparison.png"), fig; px_per_unit=3.0)
end

# Compare with IQR values for the ecoregion
begin
    ms = ["S_median" "L_median"; "S_iqr89" "L_iqr89"]
    ts = ["A) Richness" "B) Links"; "C) Richness IQR" "D) Links IQR"]
    cts = ["Expected Richness" "Expected number of links";
           "Richness 89% IQR" "Links 89% IQR"]
    cs = [:cividis :viridis; :cividis :viridis]
    fig = Figure(; resolution=(1275,600))
    for i in 1:2, j in 1:2
        m = ms[i,j]
        t = ts[i,j]
        ct = cts[i,j]
        p = background_map(fig[i,j]; title=t, titlealign=:left)
        s = surface!(ecoregion_layers["$(m)"]; colormap=cs[i,j], shading=false)
        Colorbar(p[1,2], s; height=Relative(0.5), label=ct)
    end
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_comparison_iqr.png"), fig; px_per_unit=3.0)
end

## Compare with LCBD

# Get relative LCBD values
begin
    ms = ["S" "LCBD_species"; "L" "LCBD_networks"]
    ts = ["Richness" "Species LCBD"; "Links" "Network LCBD"]
    fig = Figure(; resolution=(1275,600))
    for i in 1:2, j in 1:2
        m = ms[i,j]
        t = ts[i,j]
        p = background_map(fig[i,j]; title=t, titlealign=:left)
        s = surface!(ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        Colorbar(p[1,2], s; height=Relative(0.5), label=t)
    end
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_comparison_lcbd.png"), fig; px_per_unit=3.0)
end

# Bivariate LCBD figure for ecoregion values
function make_bivariate_figure(L1, L2, fig = Figure(); pal=bv_pal_2, kw...)
    g1 = fig[1:16, 1:4] = GridLayout()
    g2 = fig[2:5, end] = GridLayout()

    p1 = background_map(g1[1,1])
    sf = bivariatesurface!(p1, L1, L2; pal..., kw...)

    p2 = Axis(g2[1,1]; aspect = 1, xlabel = "Species LCBD", ylabel = "Network LCBD",
              xticks=0.3:0.2:0.7)
    l2 = bivariatelegend!(p2, L1, L2; pal..., kw...)
    fig
end
fig = make_bivariate_figure(
    ecoregion_layers["LCBD_species_median"],
    ecoregion_layers["LCBD_networks_median"];
    # pal=bv_pal_2,
    # rev=true,
    cmap=vec(cmat2[3:-1:1,:])
)
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_LCBD_bivariate.png"), fig; px_per_unit=3.0)
end

## Relationship between LCBD median and IQR

# Show probability densities
function make_density_figure(fig = Figure(;resolution=(800, 400)))
    ax1 = Axis(
        fig[1:3,1],
        # title="Median",
        xlabel="Relative LCBD value",
        ylabel="Probability Density"
    )
    ax2 = Axis(
        fig[1:3,2],
        # title="89% IQR",
        xlabel="89% IQR",
    )
    p1 = density!(ax1,
        unique(values(ecoregion_layers["LCBD_species_median"]));
        color=(bv_pal_2[2], 0.3),
        strokecolor=bv_pal_2[2],
        strokewidth=3,
    )
    p2 = density!(ax1,
        unique(values(ecoregion_layers["LCBD_networks_median"]));
        color=(bv_pal_2[3], 0.3),
        strokecolor=bv_pal_2[3],
        strokewidth=3,
    )
    p3 = density!(ax2,
        unique(values(ecoregion_layers["LCBD_species_iqr89"]));
        color=(bv_pal_2[2], 0.3),
        strokecolor=bv_pal_2[2],
        strokewidth=3,
    )
    p4 = density!(ax2,
        unique(values(ecoregion_layers["LCBD_networks_iqr89"]));
        color=(bv_pal_2[3], 0.3),
        strokecolor=bv_pal_2[3],
        strokewidth=3,
    )
    Legend(fig[4,:], [p1, p2], ["Species LCBD", "Network LCBD"])
    fig
end
fig = make_density_figure()
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_relation_lcbd_densities.png"), fig; px_per_unit=3.0)
end

# 3 panel version
begin
    fig = Figure(resolution=(850,1000))
    # Define layout
    g1 = fig[1:2,1] = GridLayout()
    g2 = fig[3:4,1] = GridLayout()
    g3 = fig[5,1] = GridLayout()
    # Species LCBD
    p1 = background_map(g1[1,1])
    sf1 = surface!(
        ecoregion_layers["LCBD_species_median"];
        colormap=cgrad([p0, bv_pal_2[2]]),
        shading=false
    )
    Colorbar(g1[1,2], sf1; height=Relative(0.5), label="Species LCBD")
    # Network LCBD
    p2 = background_map(g2[1,1])
    sf2 = surface!(
        ecoregion_layers["LCBD_networks_median"];
        colormap=cgrad([p0, bv_pal_2[3]]),
        shading=false
    )
    Colorbar(g2[1,2], sf2; height=Relative(0.5), label="Network LCBD")
    # Density maps
    p3 = make_density_figure(g3)
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_LCBD_all_included.png"), fig; px_per_unit=3.0)
end

# 4 panel version
begin
    fig = Figure(resolution=(1500,800))
    # Define layout
    # g1 = fig[1:4,1] = GridLayout()
    # g2 = fig[5:8,1] = GridLayout()
    # g3 = fig[1:6,2] = GridLayout()
    # g4 = fig[6:7,2] = GridLayout()
    g1 = fig[1:2,1] = GridLayout()
    g2 = fig[1:2,2] = GridLayout()
    g3 = fig[3:4,1] = GridLayout()
    g4 = fig[3:4,2] = GridLayout()
    # g3 = fig[2:4,1:3] = GridLayout()
    # g1 = fig[1:2,4:5] = GridLayout()
    # g2 = fig[3:4,4:5] = GridLayout()
    # g4 = fig[5,4:5] = GridLayout()
    # Species LCBD
    p1 = background_map(g1[1,1])
    sf1 = surface!(
        ecoregion_layers["LCBD_species_median"];
        colormap=cgrad([p0, bv_pal_2[2]]),
        shading=false
    )
    Colorbar(g1[1,2], sf1; height=Relative(0.5), label="Species LCBD")
    # Network LCBD
    p2 = background_map(g2[1,1])
    sf2 = surface!(
        ecoregion_layers["LCBD_networks_median"];
        colormap=cgrad([p0, bv_pal_2[3]]),
        shading=false
    )
    Colorbar(g2[1,2], sf2; height=Relative(0.5), label="Network LCBD")
    # Bivariate
    p3 = make_bivariate_figure(
        ecoregion_layers["LCBD_species_median"],
        ecoregion_layers["LCBD_networks_median"],
        g3;
        cmap=vec(cmat2[3:-1:1,:])
    )
    # Density maps
    p4 = make_density_figure(g4)
    # Labels
    for (label, layout) in zip(["A", "B", "C", "D"], [g1, g2, g3, g4])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            # padding = (0, 5, 5, 0),
            # halign = :right
        )
    end
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_LCBD_4panels.png"), fig; px_per_unit=3.0)
end

# Side-by-side median-median and iqr-iqr relationships
_v1 = values(ecoregion_layers["LCBD_species_median"])
_v2 = values(ecoregion_layers["LCBD_networks_median"])
_pairs_med = unique(Pair.(_v1, _v2))
_lims_med = extrema([_v1 _v2]) .+ [-0.01, 0.01]
_v3 = values(ecoregion_layers["LCBD_species_iqr89"])
_v4 = values(ecoregion_layers["LCBD_networks_iqr89"])
_pairs_iqr = unique(Pair.(_v3, _v4))
_lims_iqr = extrema([_v3 _v4]) .+ [-0.01, 0.01]
begin
    fig = Figure()
    p1 = scatter(
        fig[1,1],
        first.(_pairs_med),
        last.(_pairs_med),
        color=:black,
        axis=(;
            aspect=1,
            xlabel="Median species LCBD",
            ylabel="Median network LCBD",
            limits=(_lims_med..., _lims_med...)
        )
    )
    p2 = scatter(
        fig[1,2],
        first.(_pairs_iqr),
        last.(_pairs_iqr),
        color=:black,
        axis=(;
            aspect=1,
            xlabel="89% IQR species LCBD",
            ylabel="89% IQR network LCBD",
            limits=(_lims_iqr..., _lims_iqr...),
        )
    )
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath(fig_path, "ecoregion_relation_lcbd_iqr.png"), fig; px_per_unit=3.0)
end
