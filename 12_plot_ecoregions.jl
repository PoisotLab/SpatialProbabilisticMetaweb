#### Ecoregion plots

CAN = true
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
network_measures = ["Co", "L", "Lv", "Ld"]
measures = [network_measures..., "S", "Sv", "LCBD_species", "LCBD_networks"]
measures_ts = [
    "\nConnectance", "Number of links", "Link variance", "\nLinkage density",
    "\nRichness", "\nRichness variance", "\nRelative species LCBD", "\nRelative network LCBD"
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
    ecoregion_layers["$(o.m)_$(o.fs)"] = geotiff(
        SimpleSDMResponse, joinpath(ecoresults_path, "ecoregion_$(o.m)_$(o.fs).tif")
    )
    # Replace zero values (sites not in an ecoregion)
    ecoregion_layers["$(o.m)_$(o.fs)"] = replace(
        ecoregion_layers["$(o.m)_$(o.fs)"], 0.0 => nothing
    )
end
ecoregion_layers

## Make some plots!!

# Load worldshape shapefile to use as background on maps
ws = worldshape(50)

# Plot results
plot(
    [plot(ecoregion_layers["$(m)_median"], ws; title=m, clim=(0.0, Inf)) for m in network_measures]...;
    plot_title="Ecoregion median", xaxis="", yaxis="",
)
savefig(joinpath(fig_path, "ecoregion_all_median.png"))

# Some variations
ecoregion_plots = Dict{String, Plots.Plot}()
for (m,t) in zip(measures, measures_ts)
    ecoregion_plots[m] = plot(
        plot(ecoregion_layers["$(m)_median"], ws; c=:inferno, clim=(0.0, Inf)),
        plot(ecoregion_layers["$(m)_iqr89"], ws; c=:magma, clim=(0.0, Inf));
        layout=(2,1),
        size=(650, 600),
        cbtitle=["$t" "$t 89% IQR"],
    )
    savefig(joinpath(fig_path, "ecoregion_$m.png"))
end
ecoregion_plots["Co"]
ecoregion_plots["L"]
ecoregion_plots["Lv"]
ecoregion_plots["Ld"]
ecoregion_plots["S"]
ecoregion_plots["Sv"]
ecoregion_plots["LCBD_species"]
ecoregion_plots["LCBD_networks"]


## Compare with richness
plot(
    [plot(ecoregion_layers["$(m)_median"], ws; clim=(0.0, Inf)) for m in ["S", "Sv", "L", "Lv"]]...;
    title=["Richness" "Richness variance" "Links" "Link variance"], xaxis="", yaxis="",
)
savefig(joinpath(fig_path, "ecoregion_comparison.png"))

## Compare with LCBD

# Get relative LCBD values
plot(
    [plot(ecoregion_layers["$(m)_median"], ws; clim=(0.0, Inf)) for m in ["S", "LCBD_species", "L", "LCBD_networks"]]...;
    title=["Richness" "Species LCBD" "Links" "Network LCBD"], xaxis="", yaxis=""
)
savefig(joinpath(fig_path, "ecoregion_comparison_lcbd.png"))

## Relationship between LCBD median and IQR

# Show probability densities
p_dens = begin
    _p1 = density(unique(collect(ecoregion_layers["LCBD_species_median"])); c=bv_pal_4[3], label="Species LCBD")
    density!(unique(collect(ecoregion_layers["LCBD_networks_median"])), c=bv_pal_4[2], label="Network LCBD")
    plot!(xaxis="Relative LCBD value", yaxis="Probability density", legend=:topleft)
    _p2 = density(unique(collect(ecoregion_layers["LCBD_species_iqr89"])); c=bv_pal_4[3], label="Species LCBD")
    density!(unique(collect(ecoregion_layers["LCBD_networks_iqr89"])), c=bv_pal_4[2], label="Network LCBD")
    plot!(xaxis="89% IQR", yaxis="Probability Density")
    plot(_p1, _p2, size=(650, 400))
end
savefig(joinpath(fig_path, "ecoregion_relation_lcbd_densities.png"))

_plcbd1 = plot(ecoregion_layers["LCBD_species_median"], ws; c=cgrad([p0, bv_pal_4[3]]), cbtitle="\nRelative species LCBD");
_plcbd2 = plot(ecoregion_layers["LCBD_networks_median"], ws; c=cgrad([p0, bv_pal_4[2]]), cbtitle="\nRelative network LCBD");
begin
    _layout = @layout [a; b; [c d]]
    plot(_plcbd1, _plcbd2, _p1, _p2;
        layout=_layout, title=["a)" "b)" "c)" "d)"], size=(650, 900), titlepos=:left,
        leftmargin=2mm)
end
savefig(joinpath(fig_path, "ecoregion_LCBD_all_included.png"))

# Side-by-side median-median and iqr-iqr relationships
_v1 = collect(ecoregion_layers["LCBD_species_median"])
_v2 = collect(ecoregion_layers["LCBD_networks_median"])
_pairs_med = unique(Pair.(_v1, _v2))
_lims_med = extrema([_v1 _v2]) .+ [-0.01, 0.01]
_v3 = collect(ecoregion_layers["LCBD_species_iqr89"])
_v4 = collect(ecoregion_layers["LCBD_networks_iqr89"])
_pairs_iqr = unique(Pair.(_v3, _v4))
_lims_iqr = extrema([_v3 _v4]) .+ [-0.01, 0.01]
plot(
    scatter(first.(_pairs_med), last.(_pairs_med)),
    scatter(first.(_pairs_iqr), last.(_pairs_iqr)),
    xlab=["Median species LCBD" "89% IQR species LCBD"],
    ylab=["Median network LCBD" "89% IQR network LCBD"],
    xlims=[_lims_med _lims_iqr],
    ylims=[_lims_med _lims_iqr],
    mc=:black,
    legend=:none,
    size=(650,400),
    aspectratio=1
)
savefig(joinpath(fig_path, "ecoregion_relation_lcbd_iqr.png"))
