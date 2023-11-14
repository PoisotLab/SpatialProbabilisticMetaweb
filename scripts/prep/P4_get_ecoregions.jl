#### Ecoregions ####

CAN = true
include("../../A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    out_path = joinpath("data", "input", "canada_ecoregions.tif");
else
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    out_path = joinpath("data", "input", "quebec_ecoregions.tif")
end

# Define reference layer
reference_layer = read_geotiff(ref_path, SimpleSDMPredictor)
spatialrange = boundingbox(reference_layer)

## Ecoregions2017

# Download ecoregions from https://ecoregions.appspot.com/  (Dinerstein et al. 2017)
path = joinpath("data", "shapefiles", "ecoregions")
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
eco_file = joinpath("data", "shapefiles", "ecoregions", "Ecoregions2017.shp")
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

## Ecoregions from Canadian Government

# Download ecoregions from https://sis.agr.gc.ca/cansis/nsdb/ecostrat/gis_data.html
shp_file_can = joinpath(path, "ecoregions.shp")
if !isfile(shp_file_can)
    # Download archive
    url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip"
    zip_file = joinpath(path, "ecoregions.zip")
    isfile(zip_file) || Downloads.download(url, zip_file);
    # Unzip files
    r = ZipFile.Reader(zip_file)
    for f in r.files
        f.name == "Ecoregions/" && continue
        n = replace(f.name, "Ecoregions/" => "")
        write(joinpath(path, n), read(f, String))
    end
    # Remove archive
    rm(zip_file)
end

# Reproject to WGS84
tmp_file = string(tempname(), ".shp")
isfile(tmp_file) || run(`$(GDAL.ogr2ogr_path()) $tmp_file $shp_file_can -t_srs EPSG:4326`); # 16 min.

# Rasterize the ecoregions
eco_file_can = tmp_file
out_file_can = replace(out_path, ".tif" => "_can.tif")
sy, sx = size(reference_layer)
l, r, b, t = spatialrange
query = `$(GDAL.gdal_rasterize_path()) -a ECOREGION -init -9999 -a_nodata -9999 $eco_file_can $out_file_can -ts $sx $sy -te $l $b $r $t`;
@time run(query); # 20 sec

# Read and replace zero values
ecoregions_can = read_geotiff(out_file_can, SimpleSDMResponse)
length(unique(values(ecoregions_can))) # 194 ecoregions
extrema(ecoregions_can)

# Mask based on the Canada shapefile
ecoregions_can = mask(reference_layer, ecoregions_can)

# Investigate site difference
diff_sites = setdiff(keys(reference_layer), keys(ecoregions_can))
if length(diff_sites) > 0
    # Get site differences
    diff_layer = similar(ecoregions_can)
    diff_layer[diff_sites] = ones(length(diff_sites))

    # Map them
    CairoMakie.activate!()
    begin
        fig = heatmap(ecoregions_can)
        scatter!(diff_sites; color=:red)
        fig
    end
    save(joinpath("xtras", "ecoregioncan_mismatch_scatter.png"), fig; px_per_unit=4.0)
    fig = heatmap(diff_layer)
    save(joinpath("xtras", "ecoregioncan_mismatch_heatmap.png"), fig; px_per_unit=4.0)
    GLMakie.activate!()

    # Update reference layer
    # reference_layer = convert(SimpleSDMResponse, reference_layer)
    # reference_layer[diff_sites] = fill(nothing, length(diff_sites))
    # write_geotiff(ref_path, reference_layer) # leave as-is for now
end

# Reexport ecoregions
write_geotiff(out_file_can, convert(Float64, ecoregions_can))
