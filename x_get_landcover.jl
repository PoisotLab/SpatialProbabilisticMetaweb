include("A0_required.jl");

# Option to run for CAN
CAN = true
if (@isdefined CAN) && CAN == true
    res = 2.5;
    input_path = joinpath("data", "input");
    @info "Running for Canada at 2.5 arcmin resolution"
else
    res = 10.0;
    input_path = joinpath("xtras", "input");
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
reference_layer = read_geotiff(
    joinpath(input_path, "chelsa2_stack.tif"), SimpleSDMPredictor; bandnumber=1
)
spatialrange = boundingbox(reference_layer)

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

# Make sure to process the files in order
lc_ids = collect(eachindex(lc_files))
lc_files_ordered = [
    joinpath.(lc_path, "consensus_full_class_$id.tif") for id in lc_ids
]

# Coarsen the files
lc_layers = SimpleSDMResponse{Float32}[]
for file in lc_files_ordered
    # Coarsen to the reference resolution
    sy, sx = size(reference_layer)
    l, r, b, t = spatialrange
    out_file = tempname()
    query = `gdalwarp -r average $file $out_file -ts $sx $sy -te $l $b $r $t`
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
write_geotiff(joinpath(input_path, "landcover_stack.tif"), lc_layers)