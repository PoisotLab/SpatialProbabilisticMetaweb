include("A0_required.jl")

## Get Canada reference file

# Download shapefile
path = joinpath("data", "shapefiles", "canada")
isdir(path) || mkdir(path);
shp_file = joinpath(path, "canada.shp")
if !isfile(shp_file)
    # Download archive
    url = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000b21a_e.zip"
    zip_file = joinpath(path, "canada.zip")
    isfile(zip_file) || Downloads.download(url, zip_file);
    # Unzip files
    r = ZipFile.Reader(zip_file)
    for f in r.files
        n = replace(f.name, "lpr_000b21a_e" => "canada")
        write(joinpath(path, n), read(f, String))
    end
    # Remove archive
    rm(zip_file)
end

## Rasterize

# Reproject to WGS84
tmp_file = joinpath(path, "tmp.shp")
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
out_file = joinpath("data", "input", "canada_ref_2.tif")
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
out_file = joinpath("data", "input", "quebec_ref_10.tif")
@time run(`$(GDAL.gdal_rasterize_path()) -a PRUID $tmp_file $out_file -tr $rx $ry -te $l $b $r $t`);

# Replace zero values
test = read_geotiff(out_file, SimpleSDMResponse)
test = convert(Float32, test .∈ [10, 11, 12, 13, 24]) # Quebec & Maritimes
replace!(test, 0.0 => nothing)
write_geotiff(out_file, test)
