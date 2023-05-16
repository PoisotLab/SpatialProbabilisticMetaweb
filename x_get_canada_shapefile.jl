include("A0_required.jl")

using JSON

# From SimpleSDMLayers example
borders = download("https://raw.githubusercontent.com/AshKyd/geojson-regions/master/countries/50m/CAN.geojson")
can_data = JSON.parsefile(borders)
polys = can_data["geometry"]["coordinates"]
canada_poly = SimpleSDMLayers._format_polygon.(polys)

# Reference layer at 10 arcmin
coords = (left=-145.0, right=-50.0, bottom=40.0, top=89.0)
data_provider = RasterData(WorldClim2, BioClim)
reference_layer = SimpleSDMPredictor(data_provider; coords...)
@time reference_mask = mask(canada_poly, reference_layer) # 240 sec.
reference_mask = replace(reference_mask, 0.0 => 1.0)
heatmap(reference_mask; c=:inferno)
write_geotiff(joinpath("data", "input", "canada_ref_10.tif"), reference_mask)

# Reference layer at 2.5 arcmin
reference_layer2 = SimpleSDMPredictor(data_provider; resolution=2.5, coords...)
@time reference_mask2 = mask(canada_poly, reference_layer2) # 40 min
reference_mask2 = replace(reference_mask2, 0.0 => 1.0)
write_geotiff(joinpath("data", "input", "canada_ref_2.tif"), reference_mask2)

# Export a similar reference layer for Quebec (as in earlier analyses)
qc_coords = (left=-80.0, right=-50.0, bottom=45.0, top=65.0)
qclayer = SimpleSDMPredictor(data_provider; qc_coords...)
replace!(similar(qclayer), 0.0 => 1.0)
write_geotiff(joinpath("data", "input", "quebec_ref_10.tif"), qclayer)

## Projected canada map

#=
# Load Shapefile from Open Canada Data
canada = open("shapefiles/canada/canada.shp", "r") do io
    read(io, Shapefile.Handle)
end

# Create background map
plot(canada, c=:lightgrey, dpi=200,
     grid=false, axis=false, ticks=false,
     bg_colour=:transparent, fg_colour=:black)
savefig(joinpath("figures", "background", "x_background_ca.png"))
plot(canada, c=:lightgrey, dpi=200,
     grid=false, axis=false, ticks=false,
     bg_colour=:transparent, fg_colour=:white)
savefig(joinpath("figures", "background", "x_background_ca_white.png"))
=#