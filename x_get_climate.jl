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

## Coarsen the CHELSA layers

# Define reference layer
spatialrange = (left=-180.0, right=-40.0, bottom=18.0, top=89.0)
reference_layer = SimpleSDMPredictor(
    RasterData(WorldClim2, BioClim); resolution=res, spatialrange...
)

## Coarsen CHELSA1 layer (used to clip CHELSA2 data)

# Only one CHELSA1 layer is necessary

# Make sure we have the files locally
ch1_provider = RasterData(CHELSA1, BioClim)
_tmp = SimpleSDMPredictor(ch1_provider; layer="BIO1", spatialrange...)
_tmp = nothing

# Get the paths
ch1_path = joinpath(SimpleSDMDatasets._LAYER_PATH, "CHELSA1", "BioClim")
ch1_files = readdir(ch1_path; join=true)
filter!(contains("chelsa_bio10_01"), ch1_files)

# Coarsen the files
function coarsen_layer(file, reference_layer)
    # Coarsen to the reference resolution
    sy, sx = size(reference_layer)
    l, r, b, t = spatialrange
    out_file = tempname()
    query = `gdalwarp -r average $file $out_file -ts $sx $sy -te $l $b $r $t`
    run(query);

    # Read the coarsened file
    coarsened = read_geotiff(out_file, SimpleSDMPredictor)
end
ch1 = coarsen_layer(ch1_files[1], reference_layer)

# Export
write_geotiff(joinpath(out_path, "chelsa1_mask.tif"), ch1)

## Coarsen CHELSA2 layers

# Make sure we have the files locally
ch2_provider = RasterData(CHELSA2, BioClim)
_tmp = [
    SimpleSDMPredictor(ch2_provider; layer=l, spatialrange...) for l in layers(ch2_provider)
]
_tmp = nothing

# Get the paths
ch2_path = joinpath(SimpleSDMDatasets._LAYER_PATH, "CHELSA2", "BioClim")
ch2_files = readdir(ch2_path)

# Make sure to process the files in order
ch2_ids = collect(eachindex(ch2_files))
ch2_files_ordered = [
    joinpath.(ch2_path, "chelsa_bio$(id)_1981-2010_v.2.1.tif") for id in ch2_ids
]

# Coarsen the files
ch_layers = SimpleSDMResponse{Float32}[]
for file in ch2_files_ordered
    # Coarsen to the reference resolution
    ch2 = coarsen_layer(file, reference_layer)

    # Mask the coarsened layer with the reference layer
    ch2 = convert(Float32, ch2)
    ch2 = mask(reference_layer, ch2)

    # Warn if the sites with data differ
    diff_sites = setdiff(keys(ch2), keys(reference_layer))
    if length(diff_sites) != 0
        @warn "$(length(diff_sites)) sites that do not match between layers"
    end

    # Export the result
    push!(ch_layers, ch2)
end

# Export as a stack
write_geotiff(joinpath(out_path, "chelsa2_stack.tif"), ch_layers)
