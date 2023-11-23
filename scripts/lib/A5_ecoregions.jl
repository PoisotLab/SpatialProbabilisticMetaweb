# Load ecoregions
if (@isdefined CAN) && CAN == true
    eco_path = joinpath("data", "input", "canada_ecoregions.tif");
else
    eco_path = joinpath("data", "input", "quebec_ecoregions.tif")
end
ecoregions = read_geotiff(eco_path, SimpleSDMPredictor)

# Separate ecoregions in different layers
ecoregions_ids = unique(values(ecoregions))
ecoregions_stack = [convert(Float32, ecoregions .== id) for id in ecoregions_ids]
for e in ecoregions_stack
    replace!(e, 0.0 => nothing)
end

# Define the summary functions we will use
quantile055(x) = quantile(x, 0.055)
quantile945(x) = quantile(x, 0.945)
iqr89(x) = quantile945(x) - quantile055(x)
summary_fs = [median, iqr89]

# Define function to summarize by ecoregions
function ecoregionalize(layer, ecoregions_stack; f=median, keepzeros=true)
    l_eco = similar(layer)
    @threads for e in ecoregions_stack
        k = intersect(keys(e), keys(layer))
        l_eco[k] = fill(f(layer[k]), length(k))
    end
    if !keepzeros
        replace!(l_eco, 0.0 => nothing)
    end
    return l_eco
end