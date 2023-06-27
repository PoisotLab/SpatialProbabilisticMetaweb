#### Ecoregions ####

CAN = true
include("A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    input_path = joinpath("data", "input");
    out_path = joinpath("data", "input", "canada_ecoregions.tif");
else
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    input_path = joinpath("xtras", "input");
    out_path = joinpath("data", "input", "quebec_ecoregions.tif")
end

# Define reference layer
reference_layer = read_geotiff(ref_path, SimpleSDMPredictor)
spatialrange = boundingbox(reference_layer)

# Set the coordinates that do not match to zero
lc_layer = read_geotiff(
    joinpath(input_path, "landcover_stack.tif"), SimpleSDMPredictor; spatialrange...
)
site_mismatch = setdiff(keys(reference_layer), keys(lc_layer))
reference_layer = convert(SimpleSDMResponse, reference_layer)
reference_layer[site_mismatch] = fill(nothing, length(site_mismatch))
reference_layer = convert(SimpleSDMPredictor, reference_layer)

# Rasterize the ecoregions
eco_file = joinpath("shapefiles", "ecoregions", "Ecoregions2017.shp")
out_file = out_path
sy, sx = size(reference_layer)
l, r, b, t = spatialrange
query = `$(GDAL.gdal_rasterize_path()) -a ECO_ID $eco_file $out_file -ts $sx $sy -te $l $b $r $t`;
@time run(query); # ~ 50 sec.

# Read and replace zero values
ecoregions = read_geotiff(out_path, SimpleSDMResponse)
replace!(ecoregions, 0.0 => nothing)

# Mask based on the Canada shapefile
ecoregions = mask(reference_layer, ecoregions)

# Reexport
write_geotiff(out_path, convert(Float64, ecoregions))
