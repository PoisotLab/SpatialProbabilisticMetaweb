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

# Show relationship as scatter plot
begin
    scatter(
        unique(collect(ecoregion_layers["LCBD_species_median"])),
        unique(collect(ecoregion_layers["LCBD_species_iqr89"]));
        label="Species LCBD", c=:black
    )
    scatter!(
        unique(collect(ecoregion_layers["LCBD_networks_median"])),
        unique(collect(ecoregion_layers["LCBD_networks_iqr89"]));
        label="Network LCBD", c=:orange
    )
    plot!(xaxis=("Ecoregion median relative LCBD value", (0,1)), yaxis=("89% IQR", (0,1)))
end
savefig(joinpath(fig_path, "ecoregion_relation_lcbd_iqr.png"))

# Show probability densities
begin
    _p1 = density(unique(collect(ecoregion_layers["LCBD_species_median"])); c=:black, label="Species LCBD")
    density!(unique(collect(ecoregion_layers["LCBD_networks_median"])), c=:orange, label="Network LCBD")
    plot!(xaxis="Relative LCBD value", yaxis="Probability density")
    _p2 = density(unique(collect(ecoregion_layers["LCBD_species_iqr89"])); c=:black, label="Species LCBD")
    density!(unique(collect(ecoregion_layers["LCBD_networks_iqr89"])), c=:orange, label="Network LCBD")
    plot!(xaxis="89% IQR", yaxis="Probability Density")
    plot(_p1, _p2, size=(650, 400))
end
savefig(joinpath(fig_path, "ecoregion_relation_lcbd_densities.png"))