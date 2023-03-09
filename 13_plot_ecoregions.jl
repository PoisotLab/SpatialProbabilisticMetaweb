#### Ecoregion plots

# CAN = true
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
summary_fs = ["mean", "median", "maximum", "minimum_nonzero"]

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

# Load the ecoregion metaweb layers
ecometaweb_layers = Dict{String, SimpleSDMResponse}()
for o in opt
    ecometaweb_layers["$(o.m)_$(o.fm)"] = geotiff(
        SimpleSDMResponse, joinpath(ecoresults_path, "ecometaweb_$(o.m)_$(o.fm).tif")
    )
end

## Make some plots!!

# Plot results
plot(
    plot(ecoregion_layers["Co_mean"]; title="Co"),
    plot(ecoregion_layers["L_mean"]; title="L"),
    plot(ecoregion_layers["Lv_mean"]; title="Lv"),
    plot(ecoregion_layers["Ld_mean"]; title="Ld"),
    plot_title="Ecoregion mean"
)
savefig(joinpath(fig_path, "ecoregion_all_mean.png"))

# plot(ecoregion_layers["S"]; title="S")
# plot(ecoregion_layers["Sσ"]; title="Sσ")

# Some variations
plot(
    plot(ecoregion_layers["L_mean"]; title="mean"),
    # plot(ecoregion_layers["L_sum"]; title="sum"),
    plot(ecoregion_layers["L_median"]; title="median"),
    plot(ecoregion_layers["L_maximum"]; title="maximum"),
    plot(ecoregion_layers["L_minimum_nonzero"]; title="minimum_nonzero"),
    plot_title="L"
)
savefig(joinpath(fig_path, "ecoregion_L.png"))

## Ecoregion metaweb

# Plot results using mean of interactions
plot(
    plot(ecometaweb_layers["Co_mean"]; title="Co"),
    plot(ecometaweb_layers["L_mean"]; title="L"),
    plot(ecometaweb_layers["Lv_mean"]; title="Lv"),
    plot(ecometaweb_layers["Ld_mean"]; title="Ld"),
    plot_title="Ecoregion metaweb mean"
)
savefig(joinpath(fig_path, "ecometaweb_all_mean.png"))

# Plot results
plot(
    plot(ecometaweb_layers["L_mean"]; title="mean"),
    plot(ecometaweb_layers["L_median"]; title="median"),
    plot(ecometaweb_layers["L_maximum"]; title="maximum"),
    # plot(ecometaweb_layers["L_minimum"]; title="minimum"),
    plot(ecometaweb_layers["L_minimum_nonzero"]; title="minimum_nonzero"),
    plot_title="L ecoregion metaweb"
)
savefig(joinpath(fig_path, "ecometaweb_L.png"))