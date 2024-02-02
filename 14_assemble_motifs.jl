#### Plot motifs ####

# CAN = true
include("A0_required.jl");

# Load the corresponding results if dealing with CAN data or minimal example
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

# Make sure motif results exist
motifs_path = joinpath(results_path, "motifs")
if !isdir(motifs_path) && length(readdir(motifs_path)) < 1
    @info "No motif result. Running 13_get_motifs.jl"
    include("13_get_motifs.jl")
end

## Assemble motif results

# Check available motifs
SX_all = unique(first.(split.(readdir(motifs_path), "-")))
SX_all = replace.(SX_all, ".tif" => "")
SX_missing = setdiff(["S1", "S2", "S4", "S5"], SX_all)
if length(SX_missing) > 1
    @warn "Missing motifs $(SX_missing). Re-run 13_get_motifs.jl manually for those."
end

# Assemble target motifs
motifs = Dict{String, SimpleSDMResponse}()
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
