#### Plot motifs ####

SAVE = true
CAN = true
include("A0_required.jl");

# Load the corresponding sdm results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    results_path = joinpath("data", "results")
else
    results_path = joinpath("xtras", "results")
end

# Load motif results
motifs = Dict{String, SimpleSDMLayer}()
for SX in ["S1", "S2", "S4", "S5"]
    motifs[SX] = read_geotiff(joinpath(results_path, "$SX.tif"), SimpleSDMPredictor)
end
S1 = motifs["S1"]
S2 = motifs["S2"]
S4 = motifs["S4"]
S5 = motifs["S5"]

# Get normalized difference indexes
NDTI = (S1-S2)/(S1+S2)
NDCI = (S4-S5)/(S4+S5)

## Plot

# Plot single motifs
for (motif, layer) in motifs
    begin
        fig = background_map()
        sf = surface!(log1p(layer); shading=false)
        Colorbar(fig[1,2], sf; height=Relative(0.5), label="log($motif + 1)")
        fig
    end
    if (@isdefined SAVE) && SAVE == true
        save(joinpath("figures", "motifs_$motif.png"), fig)
    end
end

# S1-S2 comparison - Normalized difference trophic index
begin
    fig = background_map()
    sf = surface!(NDTI; shading=false, colorrange=(-1/2,1/2), colormap=:roma)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Normalized Difference Trophic Index")
    fig
end
if (@isdefined SAVE) && SAVE == true
    save(joinpath("figures", "motifs_NDI_trophic.png"), fig)
end

# S4-S5 comparison - Normalized difference competition index
begin
    fig = background_map()
    sf = surface!(NDCI; shading=false, colorrange=(-1/2,1/2), colormap=:roma)
    Colorbar(fig[1,2], sf; height=Relative(0.5), label="Normalized Difference Competition Index")
    fig
end
if (@isdefined SAVE) && SAVE == true
    save(joinpath("figures", "motifs_NDI_competition.png"), fig)
end

# Relationships
rel_pairs = ["S1" => "S2", "S4" => "S5"]
for pair in rel_pairs
    SA = motifs[pair.first]
    SB = motifs[pair.second]
    common_sites = intersect(findall(!isnothing, SA.grid), findall(!isnothing, SB.grid))
    begin
        fig = Figure()
        ax = Axis(fig[1,1],
            xlabel=("log($(pair.first)+ 1)"), ylabel=("log($(pair.second) + 1)"),
            aspect=1, yticks=0:2:10, xticks=0:2:10, limits=((nothing, 10), (nothing, 10))
        )
        hexbin!(
            log1p.(SA[common_sites]),
            log1p.(SB[common_sites]);
            bins=50
        )
        lines!(0:10, 0:10; color=:red)
        fig
    end
    if (@isdefined SAVE) && SAVE == true
        save(joinpath("figures", "motifs_rel_$(pair.first)-$(pair.second).png"), fig)
    end
end

## Ecoregion motifs figures

# Load ecoregion objects and functions
include("scripts/lib/A5_ecoregions.jl")

# Add NDI layers to motif Dict
motifs["NDTI"] = NDTI
motifs["NDCI"] = NDCI

# Create ecoregion motif layers
for SX in ["S1", "S2", "S4", "S5", "NDTI", "NDCI"], f in [mean, median, iqr89, minimum, maximum]
    motifs["$(SX)_$f"] = ecoregionalize(
        motifs["$SX"], ecoregions_stack; f=x -> f(filter(!isnan, (x)))
    )
    motifs["$(SX)_$f"] = replace(motifs["$(SX)_$f"], 0.0 => nothing)
end

# Plot ecoregion motifs
for SX in ["S1", "S2", "S4", "S5"]
    begin
        fig = background_map()
        sf = surface!(log1p(motifs["$(SX)_median"]); shading=false)
        Colorbar(fig[1,2], sf; height=Relative(0.5), label="log($SX + 1)")
        fig
    end
    if (@isdefined SAVE) && SAVE == true
        save(joinpath("figures", "ecoregions", "motifs_ecoregion_$SX.png"), fig)
    end
end

# Plot NDI variations
for (m,t) in zip(["NDTI", "NDCI"], ["trophic", "competition"]), f in [minimum, maximum]
    @info "$m $t $f"
    begin
        fig = background_map()
        sf = surface!(motifs["$(m)_$f"]; shading=false, colorrange=(-1.0, 1.0))
        Colorbar(fig[1,2], sf; height=Relative(0.5), label="Normalized Difference $(uppercasefirst(t)) Index ($f)")
        fig
    end
    if (@isdefined SAVE) && SAVE == true
        save(joinpath("figures", "ecoregions", "motifs_ecoregion_NDI_$(t)_$f.png"), fig)
    end
end
