QC = true
include("04_aggregate_sdms.jl");

# Load the previous sdm results if dealing with QC data
if (@isdefined QC) && QC == true
    fit_path = joinpath("xtras", "input", "sdm_fit_results.csv")
    results_path = joinpath("xtras", "results")
else
    fit_path = joinpath("data", "input", "sdm_fit_results.csv")
    results_path = joinpath("data", "results")
end

# Load the layers
# S = [geotiff(SimpleSDMPredictor, joinpath("data", "results", "Smeans.tif"), i) for i in 1:159]
spnames = @chain begin
    readdir(joinpath("data", "presence_absence"))
    replace.(".tif" => "")
    replace.("_" => " ")
end
D = Dict{String, SimpleSDMResponse}()
p = Progress(length(spnames))
@time @threads for i in eachindex(spnames)
    # D[spnames[i]] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "Smeans.tif"), i)
    D[spnames[i]] = geotiff(SimpleSDMPredictor, joinpath("data", "results", "Srands.tif"), i)
    # Names won't match but whatever
    # This way is faster and can run in parallel
    next!(p)
end
S = collect(values(D))

# Get LCBD values for all 4 assembly options
lcbd_layers = fill(similar(reference_layer), 4)
# for (i, S) in enumerate([Smeans, Srands, Smeans_cut, Srands_cut])
    @assert length(unique(length.(S))) == 1

    # Y matrix
    Y = reduce(hcat, collect.(S))

    # Find the negative values
    minimum(Y) # all positive for Smeans; some negative for Srands
    inds_neg = findall(<(0.0), Y) # 2 values only
    Y[inds_neg] # very very small
    minimum(Y[findall(>=(0.0), Y)]) # next smallest values are zero

    # Replace by zero
    @info "$(length(inds_neg)) negative values were replaced by zero"
    @info Y[inds_neg]
    Y[inds_neg] .= 0.0
    hellinger(Y) # now it works

    # LCBD
    lcbd_layers[i][keys(reference_layer)] = LCBD(hellinger(Y))[1]
# end

# Start over
D = nothing
S = nothing
lcbd_layers = nothing
GC.gc()