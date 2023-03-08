#### Ecoregion plots

# CAN = true
include("12_get_ecoregions_measures.jl")

## Basic network measures

# Plot results
plot(ecoregion_layers["Co"]; title="Co")
plot(ecoregion_layers["L"]; title="L")
plot(ecoregion_layers["Lv"]; title="Lv")
plot(ecoregion_layers["Ld"]; title="Ld")
# plot(ecoregion_layers["S"]; title="S")
# plot(ecoregion_layers["Sσ"]; title="Sσ")

# Some variations
plot(ecoregion_layers["S_sum"])
plot(ecoregion_layers["S_maximum"])
plot(ecoregion_layers["S_minimum"])

## Ecoregion metaweb

# Plot results using mean of interactions
plot(ecometaweb_layers["Co"]; title="Co_meta")
plot(ecometaweb_layers["L"]; title="L_meta")
plot(ecometaweb_layers["Lv"]; title="Lv_meta")
plot(ecometaweb_layers["Ld"]; title="Ld_meta")

# Plot results
plot(ecometaweb_layers["L_mean"]; title="L_meta_mean")
plot(ecometaweb_layers["L_max"]; title="L_meta_max")
plot(ecometaweb_layers["L_min"]; title="L_meta_min")