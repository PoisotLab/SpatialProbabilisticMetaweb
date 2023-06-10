#### Plot richness & LCBD results ####

CAN = true
include("A0_required.jl");

# Set corresponding resolution
if (@isdefined CAN) && CAN == true
    res = 2.5
else
    res = 10.0
end

# Load LCBD results
include("x_load_lcbd_results.jl");

# Check results
S_all
lcbd_species_all
lcbd_networks_all

## Richness plots

# Load worldshape shapefile to use as background on maps
ws = worldshape(50)

# Load background layer
spatialrange = boundingbox(S_all["mean"])
bglayer = similar(
    SimpleSDMPredictor(RasterData(WorldClim2, BioClim); resolution=2.5, spatialrange...)
)

# Richness for mean only
begin
    fig = background_map()
    hm2 = heatmap!(S_all["mean"]; colormap=:cividis)
    Colorbar(fig[1,end+1], hm2; height=Relative(0.5), label="Expected Richness")
    fig
end
save(joinpath("figures", "richness_mean.png"), fig; px_per_unit=2.0)

# GeoMakie attempt
begin
    fig = Figure()
    ga = GeoAxis(
        fig[1, 1];
        source = "+proj=longlat +datum=WGS84",
        dest = "esri:102002", # Lambert Conformal Conic
        lonlims = (spatialrange.left, spatialrange.right),
        latlims = (spatialrange.bottom, spatialrange.top),
        xlabel = "Longitude",
        ylabel = "Latitude",
    )
    hm1 = surface!(bglayer; colormap=:Greys, shading=false)
    hm2 = surface!(ga, S_all["mean"]; colormap=:cividis, shading=false)
    Colorbar(fig[1,end+1], hm2; height=Relative(0.5), label="Expected Richness")
    fig
end
save(joinpath("figures", "richness_proj_bg.png"), fig; px_per_unit=3.0)

# GeoMakie with shape
in_file = "shapefiles/ne_50m_land.shp"
out_file = "shapefiles/ne_50m_land_clip.shp"
l, r, b, t = spatialrange
query = `ogr2ogr -clipsrc $l $(b+1.0) $r $t $out_file $in_file`
run(query);
shapes = Shapefile.shapes(Shapefile.Table(out_file))
@time begin
    fig = Figure()
    ga = GeoAxis(
        fig[1, 1];
        source = "+proj=longlat +datum=WGS84",
        dest = "esri:102002", # Lambert Conformal Conic
        lonlims = (spatialrange.left, spatialrange.right),
        latlims = (spatialrange.bottom, spatialrange.top),
        xlabel = "Longitude",
        ylabel = "Latitude",
    )
    foreach(table) do geo
        poly!(ga, geo; shading=false, strokecolor=:darkgrey, strokewidth=1, color=:lightgrey)
    end
    hm2 = surface!(ga, S_all["mean"]; colormap=:cividis, shading=false)
    fig
end
save(joinpath("figures", "richness_proj_shp.png"), fig; px_per_unit=3.0)

# Richness variance for mean only
plot(Sv, ws; c=:cividis, cbtitle="Richness variance", size=(650, 400))
savefig(joinpath("figures", "richness_var.png"))

# Bivariate richness map
begin
    bivariate(S_all["mean"], Sv, ws; quantiles=true, classes=3, bv_pal_2...)
    bivariatelegend!(
        S_all["mean"],
        Sv;
        classes=3,
        # inset=(1, bbox(0.0, 0.17, 0.10, 0.28, :top, :right)),
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Expected richness",
        ylab="Richness variance",
        guidefontsize=7,
        bv_pal_2...
    )
end
savefig(joinpath("figures", "richness_bivariate.png"))

## LCBD plots

# Species LCBD
plot(lcbd_species_all["mean"], ws; c=:viridis, cbtitle="Relative species LCBD")
savefig(joinpath("figures", "lcbd_mean_species.png"))

# Network LCBD
plot(lcbd_networks_all["mean"], ws; c=:viridis, cbtitle="Relative network LCBD")
savefig(joinpath("figures", "lcbd_mean_networks.png"))

# Bivariate species-networks LCBD for mean only
begin
    bivariate(
        lcbd_networks_all["mean"], lcbd_species_all["mean"], ws;
        quantiles=true, bv_pal_4..., classes=3,
    )
    bivariatelegend!(
        lcbd_networks_all["mean"],
        lcbd_species_all["mean"];
        classes=3,
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Networks LCBD",
        ylab="Species LCBD",
        guidefontsize=7,
        bv_pal_4...
    )
end
savefig(joinpath("figures", "lcbd_bivariate_mean.png"))

# All species LCBD options
plot(
    [plot(lcbd_species_all[opt], ws; c=:viridis) for opt in options]...;
    title=titles,
    cbtitle="Relative species LCBD",
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "sampling_options", "lcbd_species_all.png"))

# All networks LCBD options
# NO DATA FOR NOW
#=
plot(
    [plot(lcbd_networks_all[opt], ws; c=:viridis) for opt in options]...;
    title=titles,
    cbtitle="Relative networks LCBD",
    layout=(2,2),
    size=(900,600),
)
savefig(joinpath("figures", "sampling_options", "lcbd_networks_all.png"))

# Bivariate species-networks LCBD
biv_plots = []
for (i, opt) in enumerate(options)
    bp = bivariate(
        lcbd_networks_all[opt], lcbd_species_all[opt], ws;
        quantiles=true, bv_pal_4..., classes=3, title=titles[i]
    )
    bp = bivariatelegend!(
        lcbd_networks_all[opt],
        lcbd_species_all[opt];
        classes=3,
        inset=(1, bbox(0.80, 0.02, 0.13, 0.28, :top, :right)),
        subplot=2,
        xlab="Networks LCBD",
        ylab="Species LCBD",
        guidefontsize=7,
        bv_pal_4...
    )
    push!(biv_plots, bp)
end
plot(biv_plots..., size = (900, 600))
savefig(joinpath("figures", "sampling_options", "lcbd_bivariate_all.png"))
=#
