using SimpleSDMLayers
using EvoTrees
using StatsBase
using StatsPlots
using ProgressMeter
using DataFrames
import CSV

# bioclim variables
layers = SimpleSDMPredictor(WorldClim, BioClim, 1:19; left=-180., right=-40., bottom=18., top=89.)
all_values = hcat([layer[keys(layer)] for layer in layers]...)

ispath(joinpath("data", "sdms")) || mkpath(joinpath("data", "sdms"))

df = [DataFrame(species = String[], occurrences = Int64[], AUC = Float64[], J = Float64[], cutoff = Float64[]) for i in 1:Threads.nthreads()]

pa_files = readdir(joinpath("data", "presence_absence"); join=true)
p = Progress(length(pa_files))

Threads.@threads for i in 1:length(pa_files)

    tree_store = EvoTreeGaussian(;
        loss=:gaussian,
        metric=:gaussian,
        nrounds=100,
        nbins=100,
        lambda=0.0,
        gamma=0.0,
        eta=0.1,
        max_depth=7,
        min_weight=1.0,
        rowsample=0.5,
        colsample=1.0,
    )

    pa_file = pa_files[i]
    spname = replace(split(pa_file, "/")[2], ".tif" => "")

    pr = geotiff(SimpleSDMResponse, pa_file, 1)
    ab = geotiff(SimpleSDMResponse, pa_file, 2)
    replace!(pr, false => nothing)
    replace!(ab, false => nothing)

    if length(pr) <= 10
        sdm = similar(layers[1])
        sdm[keys(sdm)] = fill(length(pr)/length(sdm), length(sdm))
        var = similar(sdm)
        geotiff(joinpath("data", "sdms", spname*"_model.tif"),sdm)
        geotiff(joinpath("data", "sdms", spname*"_error.tif"),var)
        next!(p)
        continue
    end

    xy_presence = keys(pr)
    xy_absence = keys(ab)
    xy = vcat(xy_presence, xy_absence)

    X = hcat([layer[xy] for layer in layers]...)
    y = vcat(fill(1.0, length(xy_presence)), fill(0.0, length(xy_absence)))
    train_size = floor(Int, 0.7 * length(y))
    train_idx = sample(1:length(y), train_size; replace=false)
    test_idx = setdiff(1:length(y), train_idx)
    Xtrain, Xtest = X[train_idx, :], X[test_idx, :]
    Ytrain, Ytest = y[train_idx], y[test_idx]
    model = fit_evotree(tree_store, Xtrain, Ytrain; X_eval=Xtest, Y_eval=Ytest)
    pred = EvoTrees.predict(model, all_values)

    distribution = similar(layers[1], Float64)
    distribution[keys(distribution)] = pred[:, 1]
    distribution

    uncertainty = similar(layers[1], Float64)
    uncertainty[keys(uncertainty)] = pred[:, 2]
    uncertainty

    cutoff = LinRange(extrema(distribution)..., 500);

    obs = y .> 0

    tp = zeros(Float64, length(cutoff));
    fp = zeros(Float64, length(cutoff));
    tn = zeros(Float64, length(cutoff));
    fn = zeros(Float64, length(cutoff));

    for (i, c) in enumerate(cutoff)
        prd = distribution[xy] .>= c
        tp[i] = sum(prd .& obs)
        tn[i] = sum(.!(prd) .& (.!obs))
        fp[i] = sum(prd .& (.!obs))
        fn[i] = sum(.!(prd) .& obs)
    end

    tpr = tp ./ (tp .+ fn);
    fpr = fp ./ (fp .+ tn);
    J = (tp ./ (tp .+ fn)) + (tn ./ (tn .+ fp)) .- 1.0;

    dx = [reverse(fpr)[i] - reverse(fpr)[i - 1] for i in 2:length(fpr)]
    dy = [reverse(tpr)[i] + reverse(tpr)[i - 1] for i in 2:length(tpr)]
    AUC = sum(dx .* (dy ./ 2.0))

    thr_index = last(findmax(J))
    τ = cutoff[thr_index]

    push!(df[Threads.threadid()], (spname, length(xy_presence), AUC, maximum(J), τ))

    range_mask = broadcast(v -> v >= τ, distribution)
    sdm = mask(range_mask, distribution)/maximum(mask(range_mask, distribution))
    geotiff(joinpath("data", "sdms", spname*"_model.tif"),distribution)
    geotiff(joinpath("data", "sdms", spname*"_error.tif"),uncertainty)
    next!(p)
end

CSV.write(joinpath("data", "input", "sdm_fit_results.csv"), vcat(df...))