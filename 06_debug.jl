# QC = true
include("A0_required.jl");


# Define paths
ref_path = joinpath("data", "input", "canada_ref_2.tif");
results_path = joinpath("data", "results")
results_options = ["Smeans", "Smeans_cut", "Srands", "Srands_cut"]
lcbd_options = ["mean", "mean_thr", "rand", "rand_thr"]

# Define common elements
spnames = @chain begin
    readdir(joinpath("data", "presence_absence"))
    replace.(".tif" => "")
    replace.("_" => " ")
end
reference_layer = geotiff(SimpleSDMPredictor, ref_path)
lcbd_layers = fill(similar(reference_layer), 4)

# Get LCBD values for all 4 assembly options
for (i, opt) in enumerate(lcbd_options[4])
    @infiltrate
    # Load the layers
    D = Dict{String, SimpleSDMResponse}()
    p = Progress(length(spnames))
    @threads for j in eachindex(spnames)
        D[spnames[j]] = geotiff(SimpleSDMPredictor, joinpath("xtras", "results", "$(results_options[4]).tif"), j)
        # Names won't match but whatever
        # This way is faster and can run in parallel
        next!(p)
    end
    S = collect(values(D))

    @assert length(unique(length.(S))) == 1

    # Y matrix
    Y = reduce(hcat, collect.(S))

    # Temporary fix for bug with negative values
    inds_neg = findall(<(0.0), Y) # 2 values only
    if length(inds_neg) > 0
        @info "$(length(inds_neg)) negative values were replaced by zero"
        @info Y[inds_neg]
        Y[inds_neg] .= 0.0
    end

    # LCBD
    lcbd_layers[i][keys(reference_layer)] = LCBD(hellinger(Y))[1]

    # Export
    geotiff(joinpath(results_path, "lcbd_species_$(opt).tif"), lcbd_layers[i])

    # Cleanup
    D = nothing
    S = nothing
    Y = nothing
    GC.gc()
end
