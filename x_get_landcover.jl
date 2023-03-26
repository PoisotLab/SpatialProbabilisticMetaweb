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
# layers = SimpleSDMPredictor(
#     WorldClim, BioClim, 1:19; resolution=res, qc_coords...
# )

# # Get the EarthEnv LandCover data
# landcover = SimpleSDMPredictor(EarthEnv, LandCover, 1:12; qc_coords...)
# landcover = [convert(Float16, l) for l in landcover]
# landcover = [replace(l, 0.0 => nothing) for l in landcover]

# # Plot 'em
# ws = worldshape(50)
# plot(layers[1], ws; title="temperature")
# plot(landcover[9], ws; title="urban")
# plot(landcover[1], ws; title="needleleaf trees")
# plot(landcover[5], ws; title="open water")
# plot(convert(Float16, broadcast(==(100.0), landcover[5])))

# # Check the values range
# extrema.(landcover)
# sort(unique(collect(landcover[1])))

## Coarsen the landcover layers

# Get the paths
lc_path = joinpath(SimpleSDMLayers._layers_assets_path, "EarthEnv", "LandCover", "partial")
lc_files = readdir(lc_path)

# Define reference layer
reference_layer = geotiff(SimpleSDMPredictor, ref_path)
spatialrange = boundingbox(reference_layer)

# Coarsen the files
lc_layers = SimpleSDMResponse{Float32}[]
for i in eachindex(lc_files)
    # Define paths (making sure 11 isn't read before 2)
    lc_file = joinpath(lc_path, "landcover_partial_$i.tif")
    out_file = tempname()

    # Coarsen to the reference resolution
    sy, sx = size(reference_layer)
    l, r, b, t = spatialrange
    query = `gdalwarp -r average $lc_file $out_file -ts $sx $sy -te $l $b $r $t`
    run(query);

    # Mask the coarsened layer with the reference layer
    lc = geotiff(SimpleSDMPredictor, out_file)
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
geotiff(joinpath(out_path, "landcover_stack.tif"), lc_layers)