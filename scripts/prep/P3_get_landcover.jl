include("../../A0_required.jl");

# Option to run for CAN
# CAN = true
if (@isdefined CAN) && CAN == true
    res = 2.5;
    input_path = joinpath("data", "input");
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    @info "Running for Canada at 2.5 arcmin resolution"
else
    res = 10.0;
    input_path = joinpath("xtras", "input");
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    @info "Running for Quebec at 10 arcmin resolution"
end

## Coarsen the landcover layers

# Define spatial range
ch2 = read_geotiff(
    joinpath(input_path, "chelsa2_stack.tif"), SimpleSDMPredictor; bandnumber=1
)
spatialrange = boundingbox(ch2)

# Load Canada reference layer for verifications later on
reference_layer = read_geotiff(ref_path, SimpleSDMPredictor)

# Get the paths
lc_path = joinpath(SimpleSDMDatasets._LAYER_PATH, "EarthEnv", "LandCover")

# Make sure we have the landcover files locally
if !isdir(lc_path)
    @info "Download EarthEnv LandCover data"
    lc_provider = RasterData(EarthEnv, LandCover)
    @threads for l in layers(lc_provider)
        SimpleSDMPredictor(lc_provider; layer=l, spatialrange...)
    end
    GC.gc()
end
lc_files = readdir(lc_path)
filter!(startswith("consensus_full_class_"), lc_files)

# Make sure to process the files in order
lc_ids = collect(eachindex(lc_files))
lc_files_ordered = [
    joinpath.(lc_path, "consensus_full_class_$id.tif") for id in lc_ids
]

# Coarsen the files
lc_layers = SimpleSDMResponse{Float32}[]
for file in lc_files_ordered
    # Coarsen to the reference resolution
    sy, sx = size(ch2)
    l, r, b, t = spatialrange
    out_file = tempname()
    query = `$(GDAL.gdalwarp_path()) -r average $file $out_file -ts $sx $sy -te $l $b $r $t`
    run(query);

    # Mask the coarsened layer with the CHELSA layer
    lc = read_geotiff(out_file, SimpleSDMPredictor)
    lc = convert(Float32, lc)
    lc = mask(ch2, lc)

    # Warn if the sites with data differ
    diff_sites = setdiff(keys(ch2), keys(lc))
    if length(diff_sites) != 0
        @warn "$(length(diff_sites)) sites that do not match between layers"
    end

    # Export the result
    push!(lc_layers, lc)
end

# Export as a stack
write_geotiff(joinpath(input_path, "landcover_stack.tif"), lc_layers)