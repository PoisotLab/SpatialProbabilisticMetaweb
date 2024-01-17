#### Ecoregion plots

SAVE = true
# CAN = true
include("A0_required.jl");

# Load the corresponding results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ecoresults_path = joinpath("data", "results", "ecoregions");
else
    ecoresults_path = joinpath("xtras", "results", "ecoregions");
end

## Basic network measures

fig_path = joinpath("figures", "ecoregions")
isdir(fig_path) || mkdir(fig_path)

# Define the network measures to use
network_measures = ["Co", "L", "Lv"]
measures = ["S", "Sv", "LCBD_species", "LCBD_networks", network_measures...]
measures_ts = [
    "Richness", "Richness variance", "Relative species LCBD", "Relative network LCBD",
    "Connectance", "Number of links", "Link variance",
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

# Put values on log scale (except for LCBD measures)
# @threads for o in opt
#     if !contains(o.m, "LCBD")
#         ecoregion_layers["$(o.m)_$(o.fs)"] = log(ecoregion_layers["$(o.m)_$(o.fs)"])
#     end
# end

# We need to fix an issue with the network LCBN layers before we compare with species LCBD
# Some sites had no links, so their LCBD values was set to nothing to avoid NaNs everywhere
# Now we'll also set them to NaN for species LCBD to compare the rest of the two layers
if length(ecoregion_layers["LCBD_species_median"]) > length(ecoregion_layers["LCBD_networks_median"])
    _nan_sites = setdiff(
        keys(ecoregion_layers["LCBD_species_median"]),
        keys(ecoregion_layers["LCBD_networks_median"])
    )
    @info "Creating a species LCBD layers without $(length(_nan_sites)) sites with missing network LCBD values"
    for f in ["median", "iqr89"]
        ecoregion_layers["LCBD_species_$f"][_nan_sites] = fill(nothing, length(_nan_sites))
    end
end

## Make some plots!!

# Single figures
ecoregion_plots = Dict{String, Figure}()
@showprogress "Ecoregion single figures:" for (m,t) in zip(measures, measures_ts)
    begin
        f = Figure(; resolution=(850,800))

        p1 = background_map(f[1,1]; title="Median", titlealign=:left)
        p2 = background_map(f[2,1]; title="89% IQR", titlealign=:left)

        sf1 = surface!(p1, ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        sf2 = surface!(p2, ecoregion_layers["$(m)_iqr89"]; colormap=:inferno, shading=false)

        Colorbar(p1[1,2], sf1; height=Relative(0.5), label="log($(t))")
        Colorbar(p2[1,2], sf2; height=Relative(0.5), label="log($(t) 89% IQR)")

        ecoregion_plots[m] = f;
    end;
    if (@isdefined SAVE) && SAVE == true
        save(joinpath(fig_path, "ecoregion_single_$m.png"), ecoregion_plots[m])
    end
end

# Check results
ecoregion_plots["S"]
ecoregion_plots["Sv"]
ecoregion_plots["LCBD_species"]
ecoregion_plots["LCBD_networks"]
ecoregion_plots["Co"]
ecoregion_plots["L"]
ecoregion_plots["Lv"]

## Compare with richness

# Compare with IQR values for the ecoregion
begin
    ms = ["S_median" "L_median"; "S_iqr89" "L_iqr89"]
    ts = ["A) Richness" "B) Links"; "C) Richness IQR" "D) Links IQR"]
    cts = ["Expected Richness" "Expected number of links";
           "Richness 89% IQR" "Links 89% IQR"]
    cm = [cgrad([p0, bv_pal_2[2]]), cgrad([p0, bv_pal_2[3]])]
    cs = ReversibleScale(log, exp)
    cticks = reshape([
        [20, 40, 60, 80], # Richness
        [10, 20, 30],     # Richness IQR
        [100, 300, 500],  # Links
        [100, 200, 300],  # Links IQR
    ], (2,2))
    fig = Figure(; resolution=(1275,600))
    for i in 1:2, j in 1:2
        m = ms[i,j]
        t = ts[i,j]
        ct = cts[i,j]
        p = background_map(fig[i,j]; title=t, titlealign=:left)
        s = surface!(ecoregion_layers["$(m)"]; colormap=cm[j], colorscale=cs, shading=false)
        Colorbar(p[1,2], s;
            height=Relative(0.5), label="$ct\n(log scale)", ticks=cticks[i,j],
            minorticksvisible=true, minorticks=IntervalsBetween(2)
        )
    end
    fig
end
if (@isdefined SAVE) && SAVE == true
    save(joinpath(fig_path, "ecoregion_comparison_iqr.png"), fig)
end

# Double bivariate for median and IQR
begin
    fig = Figure(; resolution=(800,850), figure_padding=20)

    ga = fig[1,1] = GridLayout()
    gb = fig[2,1] = GridLayout()

    # Median bivariate
    L1 = ecoregion_layers["S_median"]
    L2 = ecoregion_layers["L_median"]

    g1 = ga[1:16, 1:4] = GridLayout()
    g2 = ga[2:5, 4] = GridLayout()

    p1 = background_map(g1[1,1], title="A", titlealign=:left, titlesize=20)
    sf = bivariatesurface!(p1, L1, L2; n_stops=5, bv_pal_2...)

    p2 = Axis(g2[1,1];
        aspect = 1, xlabel = "Richness", ylabel = "Links",
        xticks=0:25:75, yticks=0:200:600
    )
    l2 = bivariatelegend!(p2, L1, L2; n_stops=5, bv_pal_2...,)

    # IQR bivariate
    L3 = ecoregion_layers["S_iqr89"]
    L4 = ecoregion_layers["L_iqr89"]

    g3 = gb[1:16, 1:4] = GridLayout()
    g4 = gb[2:5, 4] = GridLayout()

    p1 = background_map(g3[1,1], title="B", titlealign=:left, titlesize=20)
    sf = bivariatesurface!(p1, L3, L4; n_stops=5, bv_pal_2...)

    p2 = Axis(g4[1,1];
        aspect = 1, xlabel = "Richness IQR", ylabel = "Links IQR",
        xticks=0:15:45
    )
    l2 = bivariatelegend!(p2, L3, L4; n_stops=5, bv_pal_2...)

    fig
end
if (@isdefined SAVE) && SAVE == true
    save(joinpath(fig_path, "ecoregion_bivariates.png"), fig)
end

## Compare with LCBD

# Bivariate LCBD figure for ecoregion values
function make_bivariate_figure(L1, L2, fig = Figure(); pal=bv_pal_2, kw...)
    g1 = fig[1:16, 1:4] = GridLayout()
    g2 = fig[2:5, end] = GridLayout()

    p1 = background_map(g1[1,1])
    sf = bivariatesurface!(p1, L1, L2; pal..., kw...)

    p2 = Axis(g2[1,1];
        aspect = 1, xlabel = "Species LCBD", ylabel = "Network LCBD",
        xticks=0.2:0.3:0.8, yticks=0.3:0.2:0.7
    )
    l2 = bivariatelegend!(p2, L1, L2; pal..., kw...)
    fig
end
fig = make_bivariate_figure(
    ecoregion_layers["LCBD_species_median"],
    ecoregion_layers["LCBD_networks_median"];
    # pal=bv_pal_2,
    cmap=cmap2
)

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

# 4 panel version
begin
    fig = Figure(resolution=(1500,800))
    # Define layout
    g1 = fig[1:2,1] = GridLayout()
    g2 = fig[1:2,2] = GridLayout()
    g3 = fig[3:4,1] = GridLayout()
    g4 = fig[3:4,2] = GridLayout()
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
        cmap=cmap2
        # pal=bv_pal_2
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
if (@isdefined SAVE) && SAVE == true
    save(joinpath(fig_path, "ecoregion_LCBD_4panels.png"), fig)
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
if (@isdefined SAVE) && SAVE == true
    save(joinpath(fig_path, "ecoregion_relation_lcbd_iqr.png"), fig)
end
