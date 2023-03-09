#### Ecoregion plots

# CAN = true
include("12_get_ecoregions_measures.jl")

## Basic network measures

fig_path = joinpath("figures", "ecoregions")
isdir(fig_path) || mkdir(fig_path)

# Plot results
plot(
    plot(ecoregion_layers["Co"]; title="Co"),
    plot(ecoregion_layers["L"]; title="L"),
    plot(ecoregion_layers["Lv"]; title="Lv"),
    plot(ecoregion_layers["Ld"]; title="Ld"),
    plot_title="Ecoregion mean"
)
savefig(joinpath(fig_path, "ecoregion_all_mean.png"))

# plot(ecoregion_layers["S"]; title="S")
# plot(ecoregion_layers["Sσ"]; title="Sσ")

# Some variations
plot(
    plot(ecoregion_layers["L"]; title="mean"),
    # plot(ecoregion_layers["L_sum"]; title="sum"),
    plot(ecoregion_layers["L_median"]; title="median"),
    plot(ecoregion_layers["L_maximum"]; title="maximum"),
    plot(ecoregion_layers["L_minimum"]; title="minimum"),
    plot_title="L"
)
savefig(joinpath(fig_path, "ecoregion_L.png"))

## Ecoregion metaweb

# Plot results using mean of interactions
plot(
    plot(ecometaweb_layers["Co"]; title="Co"),
    plot(ecometaweb_layers["L"]; title="L"),
    plot(ecometaweb_layers["Lv"]; title="Lv"),
    plot(ecometaweb_layers["Ld"]; title="Ld"),
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