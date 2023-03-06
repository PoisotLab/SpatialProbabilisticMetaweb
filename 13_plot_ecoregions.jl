#### Ecoregion plots

# CAN = true
include("12_get_ecoregions_measures.jl")

## Basic network measures

# Plot results
plot(Co_eco; title="Co")
plot(L_eco; title="L")
plot(Lv_eco; title="Lv")
plot(Ld_eco; title="Ld")
plot(S_eco; title="S")
plot(Sσ_eco; title="Sσ")

# Some variations
plot(ecoregionalize(S, ecoregions_stack; f=sum))
plot(ecoregionalize(S, ecoregions_stack; f=maximum))
plot(ecoregionalize(S, ecoregions_stack; f=minimum))

## Ecoregion metaweb

# Plot results using mean of interactions
plot(Co_meta_eco; title="Co_meta")
plot(L_meta_eco; title="L_meta")
plot(Lv_meta_eco; title="Lv_meta")
plot(Ld_meta_eco; title="Ld_meta")

# Plot results
plot(L_meta_eco_mean; title="L_meta_mean")
plot(L_meta_eco_max; title="L_meta_max")
plot(L_meta_eco_min; title="L_meta_min")