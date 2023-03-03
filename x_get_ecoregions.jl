#### Ecoregions ####

CAN = true
include("A0_required.jl")

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    out_path = joinpath("data", "input", "canada_ecoregions.tif");
else
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    out_path = joinpath("data", "input", "quebec_ecoregions.tif")
end

# Define reference layer
reference_layer = geotiff(SimpleSDMPredictor, ref_path)
spatialrange = boundingbox(reference_layer)

# Rasterize the ecoregions
eco_file = joinpath("shapefiles", "ecoregions", "Ecoregions2017.shp")
out_file = out_path
sy, sx = size(reference_layer)
l, r, b, t = spatialrange
query = `gdal_rasterize -a ECO_ID $eco_file $out_file -ts $sx $sy -te $l $b $r $t`
@time run(query);

# Read and replace zero values
ecoregions = geotiff(SimpleSDMPredictor, out_path)
ecoregions = replace(ecoregions, 0.0 => nothing)

# Mask based on the Canada shapefile
ecoregions = mask(reference_layer, ecoregions)

# Reexport
geotiff(out_path, ecoregions)
