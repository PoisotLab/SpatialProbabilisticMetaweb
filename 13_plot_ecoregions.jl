#### Ecoregion plots

CAN = true
include("A0_required.jl")

# Load the corresponding results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ecoresults_path = joinpath("data", "ecoregions");
else
    ecoresults_path = joinpath("xtras", "ecoregions");
end

## Basic network measures

fig_path = joinpath("figures", "ecoregions")
isdir(fig_path) || mkdir(fig_path)

# Define the network measures to use
network_measures = ["Co", "L", "Lv", "Ld"]
measures = [network_measures..., "S", "Sσ"]
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
end
ecoregion_layers

## Make some plots!!

# Plot results
plot(
    [plot(ecoregion_layers["$(m)_median"]; title=m, clim=(0.0, Inf)) for m in network_measures]...,
    plot_title="Ecoregion median"
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
            [plot(ecoregion_layers["$(m)_$f"]; clim=clims) for f in summary_fs]...;
            title=permutedims([t for t in summary_ts]),
            plot_title=m
        )
    end
    savefig(joinpath(fig_path, "ecoregion_$m.png"))
end
ecoregion_plots["Co"]
ecoregion_plots["L"]
ecoregion_plots["Lv"]
ecoregion_plots["Ld"]
ecoregion_plots["S"]
ecoregion_plots["Sσ"]

## Compare with richness
plot(
    [plot(ecoregion_layers["$(m)_median"]; clim=(0.0, Inf)) for m in ["L", "Lv", "S", "Sσ"]]...;
    title=["Links" "Link variance" "Richness" "Richness variance"],
)
savefig(joinpath(fig_path, "ecoregion_comparison.png"))