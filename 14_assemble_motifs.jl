#### Plot motifs ####

# CAN = true
include("A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

## Assemble motif results

# Assemble target motifs
motifs = Dict{String, SimpleSDMResponse}()
SX_all = unique(first.(split.(readdir(joinpath("data", "results", "motifs")), "-")))
for SX in SX_all
    # Load the files for the motif
    @info "Loading motif $SX"
    files = readdir(joinpath(results_path, "motifs"); join=true)
    filter!(contains(SX), files)
    all_layers = [read_geotiff(f, SimpleSDMResponse) for f in files]

    # Make sure that all sites are unique
    all_keys = mapreduce(keys, vcat, all_layers)
    length(unique(all_keys)) == length(all_keys) || error("Non unique keys with motif $SX");

    # Create the motif layer
    SX_layer = similar(all_layers[1])
    for l in all_layers
        SX_layer[keys(l)] = values(l)
    end

    # Export the layer
    write_geotiff(joinpath(results_path, "$SX.tif"), SX_layer)

    # Empty memory
    SX_layer = nothing
    all_layers = nothing
    GC.gc()
end
GC.gc()

# Load motif results
S1 = read_geotiff(joinpath(results_path, "S1.tif"), SimpleSDMPredictor)
S2 = read_geotiff(joinpath(results_path, "S2.tif"), SimpleSDMPredictor)
S4 = read_geotiff(joinpath(results_path, "S4.tif"), SimpleSDMPredictor)
S5 = read_geotiff(joinpath(results_path, "S5.tif"), SimpleSDMPredictor)
