#### Plot motifs ####

CAN = true
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
S4 = read_geotiff(joinpath(results_path, "S4.tif"), SimpleSDMPredictor)
S5 = read_geotiff(joinpath(results_path, "S5.tif"), SimpleSDMPredictor)

## Plot

# S4
begin
    fig = background_map()
    sf = surface!(log1p(S4); shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="log(S4 + 1)")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "motifs_S4.png"), fig)
end

# S5
begin
    fig = background_map()
    sf = surface!(log1p(S5); shading=false)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="log(S5 + 1)")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "motifs_S5.png"), fig)
end

# Relationship
common_sites = intersect(findall(!isnothing, S4.grid), findall(!isnothing, S5.grid))
begin
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel=("log(S4 + 1)"), ylabel=("log(S5 + 1)"),
        aspect=1, yticks=0:2:10, xticks=0:2:10, limits=((nothing, 10), (nothing, 10))
    )
    hexbin!(
        log1p.(S4[common_sites]),
        log1p.(S5[common_sites]);
        bins=50
    )
    lines!(0:10, 0:10; color=:red)
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "motifs_S4-S5-rel.png"), fig)
end

# Log ratio
begin
    fig = background_map()
    sf = surface!(log1p(S4/S5); shading=false, colormap=:broc)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Log ratio log1p(S4/S5)")
    fig
end
if Makie.current_backend() == CairoMakie
    save(joinpath("figures", "motifs_S4-S5-logratio.png"), fig)
end

# Additional attempt
begin
    fig = background_map()
    sf = surface!(log((S4+1)/(S5+1)); shading=false, colormap=:jet)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Log ratio log((S4+1) / (S5+1))")
    fig
end

# Remove sites with zero
sites_zeros = union(findall(iszero, S4), findall(iszero, S5))
S4_nozero = convert(SimpleSDMResponse, copy(S4))
S4_nozero[sites_zeros] = fill(nothing, length(sites_zeros))
S5_nozero = convert(SimpleSDMResponse, copy(S5))
S5_nozero[sites_zeros] = fill(nothing, length(sites_zeros))
begin
    fig = background_map()
    sf = surface!(log1p(S4_nozero/S5_nozero); shading=false, colormap=:jet)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Log ratio log1p(S4/S5)")
    fig
end
