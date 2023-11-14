#### Plot motifs ####

CAN = false
include("A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

# Load CairoMakie if exporting figures
if (@isdefined SAVE) && SAVE == true
    CairoMakie.activate!()
end

# Load motif results
S3 = read_geotiff(joinpath(results_path, "S3.tif"), SimpleSDMPredictor)
S4 = read_geotiff(joinpath(results_path, "S4.tif"), SimpleSDMPredictor)


## Plot

# With background
begin
    fig = background_map(; lims=boundingbox(S3))
    sf = surface!(log(S3/S4); shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="S3")
    fig
end

# Without background
begin
    fig, ax, hm = heatmap(S3)
    Colorbar(fig[1,2], hm; height = Relative(0.5), label="S3")
    fig
end
begin
    fig, ax, hm = heatmap(S4)
    Colorbar(fig[1,2], hm; height = Relative(0.5), label="S4")
    fig
end
save(joinpath("figures", "S4.png"), fig)
begin
    fig, ax, hm = heatmap(S4; colorscale=Makie.pseudolog10)
    Colorbar(fig[1,2], hm; height = Relative(0.5), label="log(S4)")
    fig
end
begin
    fig, ax, hm = heatmap(S4; colormap=cgrad(:viridis, scale=:log10))
    Colorbar(fig[1,2], hm; height = Relative(0.5), label="log(S4)")
    fig
end
save(joinpath("figures", "S4_log.png"), fig)
begin
    fig, ax, hm = heatmap(log(S3+1)/log(S4+1))
    Colorbar(fig[1,2], hm; height = Relative(0.5), label="Log ratio")
    fig
end
begin
    fig, ax, hm = heatmap(log(S3/S4+1))
    Colorbar(fig[1,2], hm; height = Relative(0.5), label="Log ratio")
    fig
end



heatmap(S4)
heatmap(S3/S4)

heatmap(log(S3))
heatmap(log(S4))
heatmap(log(S3)/log(S4))
heatmap(log(S3/S4))