include("../../A0_required.jl")

# Make sure paths exist
shp_path = joinpath("data", "shapefiles", "canada")
isdir(shp_path) || mkdir(shp_path);

## Get Canada reference file

# Download shapefile
shp_file = joinpath(shp_path, "canada.shp")
if !isfile(shp_file)
    # Download archive
    url = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000b21a_e.zip"
    zip_file = joinpath(shp_path, "canada.zip")
    isfile(zip_file) || Downloads.download(url, zip_file);
    # Unzip files
    r = ZipFile.Reader(zip_file)
    for f in r.files
        n = replace(f.name, "lpr_000b21a_e" => "canada")
        write(joinpath(shp_path, n), read(f, String))
    end
    # Remove archive
    rm(zip_file)
end

## Rasterize

# Reproject to WGS84
tmp_file = joinpath(shp_path, "tmp.shp")
isfile(tmp_file) || run(`$(GDAL.ogr2ogr_path()) $tmp_file $shp_file -t_srs EPSG:4326`); # 16 min.

# Get file coordinates
dataset = ArchGDAL.read(tmp_file)
layer = ArchGDAL.getlayer(dataset, 0)
envelope = ArchGDAL.envelope(layer)
rx, ry = (2.5/60.0, 2.5/60.0) # 2.5 arcmin
# coords = (
#     left = floor(envelope.MinX/rx)*rx,
#     right = ceil(envelope.MaxX/rx)*rx,
#     bottom = floor(envelope.MinY/rx)*rx,
#     top = ceil(envelope.MaxY/rx)*rx,
# )
coords = (
    left = floor(envelope.MinX),
    right = ceil(envelope.MaxX),
    bottom = floor(envelope.MinY),
    top = ceil(envelope.MaxY),
)
l, r, b, t = coords

# Rasterize
out_path = joinpath("data", "input")
ispath(out_path) || mkpath(out_path)
out_file = joinpath(out_path, "canada_ref_2.tif")
@time run(`$(GDAL.gdal_rasterize_path()) -a PRUID $tmp_file $out_file -tr $rx $ry -te $l $b $r $t`);

# Replace zero values
test = read_geotiff(out_file, SimpleSDMResponse)
replace!(test, 0.0 => nothing)
test = convert(Float32, test .> 0.0)
write_geotiff(out_file, test)

## Repeat with smaller layer for Quebec & Maritimes

# Get file coordinates
rx, ry = (10/60, 10/60) # 10 arcmin
l, r, b, t = (left=-80.0, right=-50.0, bottom=45.0, top=65.0)

# Rasterize
out_file = joinpath(out_path, "quebec_ref_10.tif")
@time run(`$(GDAL.gdal_rasterize_path()) -a PRUID $tmp_file $out_file -tr $rx $ry -te $l $b $r $t`);

# Replace zero values
test = read_geotiff(out_file, SimpleSDMResponse)
test = convert(Float32, test .∈ [10, 11, 12, 13, 24]) # Quebec & Maritimes
replace!(test, 0.0 => nothing)
write_geotiff(out_file, test)

## Even smaller layer for New Brunswick

# Get file coordinates
rx, ry = (10/60, 10/60) # 10 arcmin
l, r, b, t = (left=-80.0, right=-50.0, bottom=45.0, top=65.0)

# Rasterize
out_file = joinpath(out_path, "newbrunswick_ref_10.tif")
@time run(`$(GDAL.gdal_rasterize_path()) -a PRUID $tmp_file $out_file -tr $rx $ry -te $l $b $r $t`);

# Replace zero values
test = read_geotiff(out_file, SimpleSDMResponse)
test = convert(Float32, test .∈ [13]) # Quebec & Maritimes
replace!(test, 0.0 => nothing)
test = clip(test; left=-69.0, right=-63.0, top=49.0)
write_geotiff(out_file, test)