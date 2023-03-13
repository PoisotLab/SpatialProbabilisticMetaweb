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
network_fs = [connectance, links, links_var, linkage_density]
network_filename = ["connectance", "links_mean", "links_var", "links_density"]
summary_fs = ["median", "quantile055", "quantile945", "iqr89"]
summary_ts = ["median", "5.5% quantile", "94.5% quantile", "89% IQR"]

# Predefine set of options
opt = []
for (m, fn) in zip(network_measures, network_fs), fm in summary_fs
    o = (m = m, fn = fn, fm = fm)
    push!(opt, o)
end

# Load the ecoregion summary layers
ecoregion_layers = Dict{String, SimpleSDMResponse}()
for o in opt
    ecoregion_layers["$(o.m)_$(o.fm)"] = geotiff(
        SimpleSDMResponse, joinpath(ecoresults_path, "ecoregion_$(o.m)_$(o.fm).tif")
    )
end

## Make some plots!!

# Plot results
plot(
    [plot(ecoregion_layers["$(m)_median"]; title=m) for m in network_measures]...,
    plot_title="Ecoregion median"
)
savefig(joinpath(fig_path, "ecoregion_all_median.png"))

# plot(ecoregion_layers["S"]; title="S")
# plot(ecoregion_layers["Sσ"]; title="Sσ")

# Some variations
ecoregion_plots = Dict{String, Plots.Plot}()
for m in network_measures
    begin
        _L = [ecoregion_layers["$(m)_$f"] for f in summary_fs]
        clim1 = mapreduce(minimum, min, _L)
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
