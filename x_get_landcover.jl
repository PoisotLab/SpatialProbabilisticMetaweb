include("A0_required.jl");

# Option to run for CAN
CAN = true
if (@isdefined CAN) && CAN == true
    res = 2.5;
    ref_path = joinpath("data", "input", "canada_ref_2.tif");
    out_path = joinpath("data", "input");
    @info "Running for Canada at 2.5 arcmin resolution"
else
    res = 10.0;
    ref_path = joinpath("data", "input", "quebec_ref_10.tif");
    out_path = joinpath("xtras", "input");
    @info "Running for Quebec at 10 arcmin resolution"
end

# # Load all BIOCLIM variables
# qc_coords = (left=-80.0, right=-50.0, bottom=45.0, top=65.0)
# wc_provider = RasterData(WorldClim2, BioClim)
# wc_layers = [
#     SimpleSDMPredictor(wc_provider; layer=l, resolution=res, qc_coords...) for
#     l in layers(wc_provider)
# ]

# # Get the EarthEnv LandCover data
# lc_provider = RasterData(EarthEnv, LandCover)
# lc_layers = [
#     SimpleSDMPredictor(lc_provider; layer=l, qc_coords...) for l in layers(lc_provider)
# ]
# lc_layers = [convert(Float16, l) for l in lc_layers]
# lc_layers = [replace(l, 0.0 => nothing) for l in lc_layers]

# # Plot 'em
# heatmap(wc_layers[1]; colormap=:inferno, axis=(; title="temperature"))
# heatmap(lc_layers[9]; colormap=:inferno, axis=(; title="urban"))
# heatmap(lc_layers[1]; colormap=:inferno, axis=(; title="needleleaf trees"))
# heatmap(lc_layers[5]; colormap=:inferno, axis=(; title="open water"))
# heatmap(lc_layers[5] .== 100; colormap=:inferno, axis=(; title="open water"))

# # Check the values range
# extrema.(lc_layers)
# sort(unique(values(lc_layers[1])))

## Coarsen the landcover layers

# Define reference layer
spatialrange = (left=-180.0, right=-40.0, bottom=18.0, top=89.0)
reference_layer = SimpleSDMPredictor(
    RasterData(WorldClim2, BioClim); resolution=res, spatialrange...
)

# Make sure we have the landcover files locally
lc_provider = RasterData(EarthEnv, LandCover)
_tmp = [
    SimpleSDMPredictor(lc_provider; layer=l, spatialrange...) for l in layers(lc_provider)
]
_tmp = nothing

# Get the paths
lc_path = joinpath(SimpleSDMDatasets._LAYER_PATH, "EarthEnv", "LandCover")
lc_files = readdir(lc_path)
filter!(startswith("consensus_full_class_"), lc_files)

# Coarsen the files
lc_ids = replace.(lc_files, "consensus_full_class_" => "", ".tif" => "")
lc_layers = SimpleSDMResponse{Float32}[]
for id in lc_ids
    # Define paths (making sure 11 isn't read before 2)
    lc_file = joinpath(lc_path, "consensus_full_class_$id.tif")
    out_file = tempname()

    # Coarsen to the reference resolution
    sy, sx = size(reference_layer)
    l, r, b, t = spatialrange
    query = `gdalwarp -r average $lc_file $out_file -ts $sx $sy -te $l $b $r $t`
    run(query);

    # Mask the coarsened layer with the reference layer
    lc = read_geotiff(out_file, SimpleSDMPredictor)
    lc = convert(Float32, lc)
    lc = mask(reference_layer, lc)

    # Warn if the sites with data differ
    diff_sites = setdiff(keys(lc), keys(reference_layer))
    if length(diff_sites) != 0
        @warn "$(length(diff_sites)) sites that do not match between layers"
    end

    # Export the result
    push!(lc_layers, lc)
end

# Export as a stack
write_geotiff(joinpath(out_path, "landcover_stack.tif"), lc_layers)