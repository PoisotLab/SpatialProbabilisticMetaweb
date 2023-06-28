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

## Ecoregions2017

# Download ecoregions from https://ecoregions.appspot.com/  (Dinerstein et al. 2017)
path = joinpath("shapefiles", "ecoregions")
isdir(path) || mkdir(path);
shp_file = joinpath(path, "Ecoregions2017.shp")
if !isfile(shp_file)
    # Download archive
    url = "https://storage.googleapis.com/teow2016/Ecoregions2017.zip"
    zip_file = joinpath(path, "Ecoregions2017.zip")
    isfile(zip_file) || Downloads.download(url, zip_file);
    # Unzip files
    r = ZipFile.Reader(zip_file)
    for f in r.files
        write(joinpath(path, f.name), read(f, String))
    end
    # Remove archive
    rm(zip_file)
end

# Rasterize the ecoregions
eco_file = joinpath("shapefiles", "ecoregions", "Ecoregions2017.shp")
out_file = out_path
sy, sx = size(reference_layer)
l, r, b, t = spatialrange
query = `$(GDAL.gdal_rasterize_path()) -a ECO_ID -init -9999 -a_nodata -9999 $eco_file $out_file -ts $sx $sy -te $l $b $r $t`;
@time run(query); # ~ 50 sec.

# Read and replace zero values
ecoregions = read_geotiff(out_file, SimpleSDMResponse)
length(unique(values(ecoregions))) # 45 ecoregions
extrema(ecoregions) # CAREFUL: One ecoregion ID is 0 (and it's normal)

# Mask based on the Canada shapefile
ecoregions = mask(reference_layer, ecoregions)

# Investigate site difference
diff_sites = setdiff(keys(reference_layer), keys(ecoregions))
if length(diff_sites) > 0
    # Get site differences
    diff_layer = similar(ecoregions)
    diff_layer[diff_sites] = ones(length(diff_sites))

    # Map them
    CairoMakie.activate!()
    begin
        fig = heatmap(ecoregions)
        scatter!(diff_sites; color=:red)
        fig
    end
    save(joinpath("xtras", "ecoregion_mismatch_scatter.png"), fig; px_per_unit=4.0)
    fig = heatmap(diff_layer)
    save(joinpath("xtras", "ecoregion_mismatch_heatmap.png"), fig; px_per_unit=4.0)
    GLMakie.activate!()

    # Update reference layer
    # reference_layer = convert(SimpleSDMResponse, reference_layer)
    # reference_layer[diff_sites] = fill(nothing, length(diff_sites))
    # write_geotiff(ref_path, reference_layer) # leave as-is for now
end

# Reexport ecoregions
write_geotiff(out_file, convert(Float64, ecoregions))

