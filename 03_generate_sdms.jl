include("A0_required.jl")

# Option to run for CAN
# CAN = true
if (@isdefined CAN) && CAN == true
    res = 2.5;
    pa_path = joinpath("data", "presence_absence");
    sdm_path = joinpath("data", "sdms");
    input_path = joinpath("data", "input");
    @info "Running for Canada at 2.5 arcmin resolution"
else
    res = 10.0;
    pa_path = joinpath("xtras", "presence_absence");
    sdm_path = joinpath("xtras", "sdms");
    input_path = joinpath("xtras", "input");
    @info "Running for Quebec at 10 arcmin resolution"
end

# Load all BIOCLIM and EarthEnv variables
wc_layers = SimpleSDMPredictor(
    WorldClim, BioClim, 1:19; resolution=res, left=-180.0, right=-40.0, bottom=18.0, top=89.0
)
lc_layers = [
    geotiff(SimpleSDMPredictor, joinpath(input_path, "landcover_stack.tif"), i) for i in 1:12
]

# Set the coordinates that do not match to zero
site_mismatch = setdiff(keys(wc_layers[1]), keys(lc_layers[1]))
wc_layers = convert.(SimpleSDMResponse, wc_layers)
for wc in wc_layers
    wc[site_mismatch] = fill(nothing, length(site_mismatch))
end
wc_layers = convert.(SimpleSDMPredictor, wc_layers)

# Assemble all layers
layers = [wc_layers..., lc_layers...]
all_values = mapreduce(collect, hcat, layers)

# Verify that output path exists
isdir(sdm_path) || mkpath(sdm_path)

# Empty DataFrame to collect model statistics
df = [DataFrame(species = String[], occurrences = Int64[], AUC = Float64[], J = Float64[], cutoff = Float64[]) for i in 1:Threads.nthreads()]

# List all species files
pa_files = readdir(pa_path; join=true)
filter!(contains(".tif"), pa_files)

# Run SDMs, one species per loop
p = Progress(length(pa_files))
Threads.@threads for i in 1:length(pa_files)
    # Seed for reproducibility
    Random.seed!(i)

    # Model parameters
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

    # Species name
    pa_file = pa_files[i]
    spname = replace(split(pa_file, "/")[3], ".tif" => "")

    # Presence-absence files
    pr = geotiff(SimpleSDMResponse, pa_file, 1)
    ab = geotiff(SimpleSDMResponse, pa_file, 2)
    replace!(pr, false => nothing)
    replace!(ab, false => nothing)

    # Exit if too few presence sites
    if length(pr) <= 10
        sdm = similar(layers[1])
        sdm[keys(sdm)] = fill(length(pr)/length(sdm), length(sdm))
        var = similar(sdm)
        geotiff(joinpath(sdm_path, spname*"_model.tif"),sdm)
        geotiff(joinpath(sdm_path, spname*"_error.tif"),var)
        if !(@isdefined quiet) || quiet == false
            # Print progress bar
            next!(p)
        end
        continue
    end

    # Coordinates for presence & absence sites
    xy_presence = keys(pr)
    xy_absence = keys(ab)
    xy = vcat(xy_presence, xy_absence)
    xy_inds = indexin(xy, keys(layers[1]))

    # Assemble data for models
    X = @view all_values[xy_inds, :]
    y = vcat(fill(1.0, length(xy_presence)), fill(0.0, length(xy_absence)))
    train_size = floor(Int, 0.7 * length(y))
    train_idx = sample(1:length(y), train_size; replace=false)
    test_idx = setdiff(1:length(y), train_idx)
    Xtrain, Xtest = view(X, train_idx, :), view(X, test_idx, :)
    Ytrain, Ytest = view(y, train_idx), view(y, test_idx)

    # Make sure landcover variables are not all zero
    _all_zero = map(eachcol(X)) do col
        all(iszero.(col))
    end
    if any(isone.(_all_zero))
        _inds_all_zero = findall(isone, _all_zero)
        "$spname: columns $(_inds_all_zero) only contain zeros"
    end

    # Fit & run model
    model = fit_evotree(tree_store, Xtrain, Ytrain; X_eval=Xtest, Y_eval=Ytest)
    pred = EvoTrees.predict(model, all_values)

    # Assemble distribution layer
    distribution = similar(layers[1], Float64)
    distribution[keys(distribution)] = pred[:, 1]
    distribution

    # Assemble uncertainty layer
    uncertainty = similar(layers[1], Float64)
    uncertainty[keys(uncertainty)] = pred[:, 2]
    uncertainty

    # Summary statistics
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

    # Finalize layers
    range_mask = broadcast(v -> v >= τ, distribution)
    sdm = mask(range_mask, distribution)/maximum(mask(range_mask, distribution))
    geotiff(joinpath(sdm_path, spname*"_model.tif"),distribution)
    geotiff(joinpath(sdm_path, spname*"_error.tif"),uncertainty)

    if !(@isdefined quiet) || quiet == false
        # Print progress bar
        next!(p)
    end
end

# Export model statistics
CSV.write(joinpath(input_path, "sdm_fit_results.csv"), vcat(df...))