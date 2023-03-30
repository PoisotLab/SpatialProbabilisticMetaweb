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
summary_fs = ["median", "quantile055", "quantile945", "iqr89"]
summary_ts = ["median", "5.5% quantile", "94.5% quantile", "89% IQR"]

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
    # plot_title="Ecoregion median", size=(700,400), left_margin=3mm
    plot_title="Ecoregion median", xaxis="", yaxis="",
)
savefig(joinpath(fig_path, "ecoregion_all_median.png"))

# Some variations
ecoregion_plots = Dict{String, Plots.Plot}()
for m in measures
    begin
        _L = [ecoregion_layers["$(m)_$f"] for f in summary_fs]
        # clim1 = mapreduce(minimum, min, _L)
        clim1 = 0.0
        clim2 = mapreduce(maximum, max, _L)
        clims = (clim1, clim2)
        ecoregion_plots[m] = plot(
            [plot(ecoregion_layers["$(m)_$f"], ws; clim=clims) for f in summary_fs]...;
            title=permutedims([t for t in summary_ts]),
            plot_title=m,
            xaxis="",
            yaxis="",
        )
    end
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