include("../../A0_required.jl");

# Option to run for CAN
# CAN = true
if (@isdefined CAN) && CAN == true
    res = 2.5;
    out_path = joinpath("data", "input");
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    @info "Running for Canada at 2.5 arcmin resolution"
else
    res = 10.0;
    out_path = joinpath("xtras", "input");
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    @info "Running for Quebec at 10 arcmin resolution"
end

# Make sure paths exist
ispath(out_path) || mkpath(out_path)

## Coarsen the CHELSA layers

# Define spatial range
spatialrange = (left=-180.0, right=-40.0, bottom=18.0, top=89.0)

# Load Canada reference layer for verifications later on
reference_layer = read_geotiff(ref_path, SimpleSDMPredictor)

## Coarsen CHELSA1 layer (used to clip CHELSA2 data)

# Only one CHELSA1 layer is necessary

# Get the paths
ch1_path = joinpath(SimpleSDMDatasets._LAYER_PATH, "CHELSA1", "BioClim")

# Make sure we have the files locally
if !isdir(ch1_path)
    ch1_provider = RasterData(CHELSA1, BioClim)
    SimpleSDMPredictor(ch1_provider; layer="BIO1", spatialrange...)
end
ch1_files = readdir(ch1_path; join=true)
filter!(contains("chelsa_bio10_01"), ch1_files)

# Coarsen the files
function coarsen_layer(file, coords, res)
    # Coarsen to the reference resolution
    l, r, b, t = spatialrange
    ry, rx = res
    out_file = tempname()
    query = `$(GDAL.gdalwarp_path()) -r average $file $out_file -tr $rx $ry -te $l $b $r $t`
    run(query);

    # Read the coarsened file
    coarsened = read_geotiff(out_file, SimpleSDMPredictor)
end
ch1 = coarsen_layer(ch1_files[1], spatialrange, (res/60, res/60))

# Convert to Boolean layer
ch1 = broadcast(!isnothing, ch1)

# Export
write_geotiff(joinpath(out_path, "chelsa1_mask.tif"), convert(Float32, ch1))

## Coarsen CHELSA2 layers

# Get the paths
ch2_path = joinpath(SimpleSDMDatasets._LAYER_PATH, "CHELSA2", "BioClim")

# Make sure we have the files locally
if !isdir(ch2_path)
    ch2_provider = RasterData(CHELSA2, BioClim)
    @threads for l in layers(ch2_provider)
        SimpleSDMPredictor(ch2_provider; layer=l, spatialrange...)
    end
    GC.gc()
end
ch2_files = readdir(ch2_path)

# Make sure to process the files in order
ch2_ids = collect(eachindex(ch2_files))
ch2_files_ordered = [
    joinpath.(ch2_path, "chelsa_bio$(id)_1981-2010_v.2.1.tif") for id in ch2_ids
]

# Coarsen the files
ch_layers = SimpleSDMResponse{Float32}[]
ch_diff_sites = []
for file in ch2_files_ordered
    # Coarsen to the reference resolution
    ch2 = coarsen_layer(file, spatialrange, (res/60, res/60))

    # Mask the coarsened layer with the reference layer
    ch2 = convert(Float32, ch2)
    ch2 = mask(ch1, ch2)

    # Warn if the sites with data differ
    diff_sites = setdiff(keys(reference_layer), keys(ch2))
    if length(diff_sites) != 0
        @warn "$(length(diff_sites)) missing compared to reference_layer"
    end

    # Export the result
    push!(ch_layers, ch2)
    push!(ch_diff_sites, diff_sites)
end

# Investigate missing sites
diff_sites = unique(reduce(vcat, ch_diff_sites))
if length(diff_sites) != 0
    # Map missing sites
    begin
        fig = heatmap(ch_layers[1])
        scatter!(diff_sites; color=:red)
        fig
    end # Only few mismatches near the coasts

    # Update reference layer
    reference_layer = convert(SimpleSDMResponse, reference_layer)
    reference_layer[diff_sites] = fill(nothing, length(diff_sites))
    write_geotiff(ref_path, reference_layer)
end

# Export CHELSA layers as a stack
write_geotiff(joinpath(out_path, "chelsa2_stack.tif"), ch_layers)
