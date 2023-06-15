#### Ecoregion plots

CAN = true
include("A0_required.jl");

# Load the corresponding results if dealing with QC or CAN data
if (@isdefined CAN) && CAN == true
    ecoresults_path = joinpath("data", "results", "ecoregions");
else
    ecoresults_path = joinpath("xtras", "results", "ecoregions");
end

## Basic network measures

fig_path = joinpath("figures", "ecoregions")
isdir(fig_path) || mkdir(fig_path)

# Define the network measures to use
network_measures = ["Co", "L", "Lv", "Ld"]
measures = [network_measures..., "S", "Sv", "LCBD_species", "LCBD_networks"]
measures_ts = [
    "Connectance", "Number of links", "Link variance", "Linkage density",
    "Richness", "Richness variance", "Relative species LCBD", "Relative network LCBD"
]
summary_fs = ["median", "iqr89"]
summary_ts = ["median", "89% IQR"]

# Predefine set of options
opt = []
for m in measures, fs in summary_fs
    o = (m = m, fs = fs)
    push!(opt, o)
end
opt

# Load the ecoregion summary layers
ecoregion_layers = Dict{String, SimpleSDMResponse}()
for o in opt
    # Load layer
    path = joinpath(ecoresults_path, "ecoregion_$(o.m)_$(o.fs).tif")
    ecoregion_layers["$(o.m)_$(o.fs)"] = read_geotiff(path, SimpleSDMResponse)
    # Replace zero values (sites not in an ecoregion)
    replace!(ecoregion_layers["$(o.m)_$(o.fs)"], 0.0 => nothing)
end
ecoregion_layers

## Make some plots!!

# Plot results
begin
    fig = Figure(; resolution=(1200,600))
    for i in 1:2, j in 1:2
        m = reshape(network_measures, (2,2))[i,j]
        t = reshape(measures_ts[1:4], (2,2))[i,j]
        background_map(; fig=fig, pos=[i,j], title=t, titlealign=:left)
        hm2 = surface!(ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        Colorbar(fig[i, j][1,2], hm2; height=Relative(0.5))
    end
    fig
end
save(joinpath(fig_path, "ecoregion_all_median.png"), fig; px_per_unit=3.0)

# Some variations
ecoregion_plots = Dict{String, Figure}()
for (m,t) in zip(measures, measures_ts)
    begin
        fig = Figure(; resolution=(800,800))

        background_map(; fig=fig, pos=[1,1], title="Median", titlealign=:left)
        background_map(; fig=fig, pos=[2,1], title="89% IQR", titlealign=:left)

        hm1 = surface!(fig[1,1], ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        hm2 = surface!(fig[2,1], ecoregion_layers["$(m)_iqr89"]; colormap=:inferno, shading=false)

        Colorbar(fig[1, 1][1,2], hm1; height=Relative(0.5), label="$(t)")
        Colorbar(fig[2, 1][1,2], hm2; height=Relative(0.5), label="$(t) 89% IQR")

        ecoregion_plots[m] = fig;
    end;
end
ecoregion_plots["Co"]
ecoregion_plots["L"]
ecoregion_plots["Lv"]
ecoregion_plots["Ld"]
ecoregion_plots["S"]
ecoregion_plots["Sv"]
ecoregion_plots["LCBD_species"]
ecoregion_plots["LCBD_networks"]

# Export
@time @threads for m in String.(keys(ecoregion_plots))
    @time save(joinpath(fig_path, "ecoregion_$m.png"), ecoregion_plots[m]; px_per_unit=3.0)
end

## Compare with richness
begin
    ms = ["S" "Sv"; "L" "Lv"]
    ts = ["Richness" "Richness variance"; "Links" "Link variance"]
    fig = Figure(; resolution=(1200,600))
    for i in 1:2, j in 1:2
        m = ms[i,j]
        t = ts[i,j]
        background_map(; fig=fig, pos=[i,j], title=t, titlealign=:left)
        hm2 = surface!(ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        Colorbar(fig[i, j][1,2], hm2; height=Relative(0.5), label=t)
    end
    fig
end
save(joinpath(fig_path, "ecoregion_comparison.png"), fig; px_per_unit=3.0)

## Compare with LCBD

# Get relative LCBD values
begin
    ms = ["S" "LCBD_species"; "L" "LCBD_networks"]
    ts = ["Richness" "Species LCBD"; "Links" "Network LCBD"]
    fig = Figure(; resolution=(1200,600))
    for i in 1:2, j in 1:2
        m = ms[i,j]
        t = ts[i,j]
        background_map(; fig=fig, pos=[i,j], title=t, titlealign=:left)
        hm2 = surface!(ecoregion_layers["$(m)_median"]; colormap=:inferno, shading=false)
        Colorbar(fig[i, j][1,2], hm2; height=Relative(0.5), label=t)
    end
    fig
end
save(joinpath(fig_path, "ecoregion_comparison_lcbd.png"), fig; px_per_unit=3.0)

## Relationship between LCBD median and IQR

# Show probability densities
p_dens = begin
    _p1 = density(unique(collect(ecoregion_layers["LCBD_species_median"])); c=bv_pal_4[3], label="Species LCBD")
    density!(unique(collect(ecoregion_layers["LCBD_networks_median"])), c=bv_pal_4[2], label="Network LCBD")
    plot!(xaxis="Relative LCBD value", yaxis="Probability density", legend=:topleft)
    _p2 = density(unique(collect(ecoregion_layers["LCBD_species_iqr89"])); c=bv_pal_4[3], label="Species LCBD")
    density!(unique(collect(ecoregion_layers["LCBD_networks_iqr89"])), c=bv_pal_4[2], label="Network LCBD")
    plot!(xaxis="89% IQR", yaxis="Probability Density")
    plot(_p1, _p2, size=(650, 400))
end
savefig(joinpath(fig_path, "ecoregion_relation_lcbd_densities.png"))

_plcbd1 = plot(ecoregion_layers["LCBD_species_median"], ws; c=cgrad([p0, bv_pal_4[3]]), cbtitle="\nRelative species LCBD");
_plcbd2 = plot(ecoregion_layers["LCBD_networks_median"], ws; c=cgrad([p0, bv_pal_4[2]]), cbtitle="\nRelative network LCBD");
begin
    _layout = @layout [a; b; [c d]]
    plot(_plcbd1, _plcbd2, _p1, _p2;
        layout=_layout, title=["a)" "b)" "c)" "d)"], size=(650, 900), titlepos=:left,
        leftmargin=2mm)
end
savefig(joinpath(fig_path, "ecoregion_LCBD_all_included.png"))

# Side-by-side median-median and iqr-iqr relationships
_v1 = values(ecoregion_layers["LCBD_species_median"])
_v2 = values(ecoregion_layers["LCBD_networks_median"])
_pairs_med = unique(Pair.(_v1, _v2))
_lims_med = extrema([_v1 _v2]) .+ [-0.01, 0.01]
_v3 = values(ecoregion_layers["LCBD_species_iqr89"])
_v4 = values(ecoregion_layers["LCBD_networks_iqr89"])
_pairs_iqr = unique(Pair.(_v3, _v4))
_lims_iqr = extrema([_v3 _v4]) .+ [-0.01, 0.01]
plot(
    scatter(first.(_pairs_med), last.(_pairs_med)),
    scatter(first.(_pairs_iqr), last.(_pairs_iqr)),
    xlab=["Median species LCBD" "89% IQR species LCBD"],
    ylab=["Median network LCBD" "89% IQR network LCBD"],
    xlims=[_lims_med _lims_iqr],
    ylims=[_lims_med _lims_iqr],
    mc=:black,
    legend=:none,
    size=(650,400),
    aspectratio=1
)
savefig(joinpath(fig_path, "ecoregion_relation_lcbd_iqr.png"))
