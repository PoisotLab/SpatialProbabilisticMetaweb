include("A0_required.jl")

# Load Shapefile from Open Canada Data
canada = open("shapefiles/canada/lpr_000b16a_e.shp", "r") do io
    read(io, Shapefile.Handle)
end

# Create background map
plot(canada, c=:lightgrey, dpi=200,
     grid=false, axis=false, ticks=false,
     bg_colour=:transparent, fg_colour=:black)
savefig("figures", "x_background_ca.png")
plot(canada, c=:lightgrey, dpi=200,
     grid=false, axis=false, ticks=false,
     bg_colour=:transparent, fg_colour=:white)
savefig("figures", "x_background_ca_white.png")


# Alternative
# From SimpleSDMLayers example
borders = download("https://raw.githubusercontent.com/AshKyd/geojson-regions/master/countries/50m/CAN.geojson")
can_data = JSON.parsefile(borders)
polys = can_data["geometry"]["coordinates"]
canada_poly = SimpleSDMLayers._format_polygon.(polys)

# Reference layer at 10 arcmin
coords = (left=-145.0, right=-50.0, bottom=40.0)
reflayer = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
@time reflayer_mask = mask(canada_poly, reflayer) # 140 sec.
plot(reflayer_mask)
geotiff(joinpath("data", "input", "canada_ref_10.tif"), similar(reflayer_mask))

# Reference layer at 2.5 arcmin
reflayer2 = SimpleSDMPredictor(WorldClim, BioClim, 1; resolution=2.5, coords...)
@time reflayer2_mask = mask(canada_poly, reflayer2) # 40 min
geotiff(joinpath("data", "input", "canada_ref_2.tif"), similar(reflayer2_mask))

## Testing on the model layers
# Clip the layers to Canada
layer = geotiff(SimpleSDMPredictor, "data/sdms/Aeorestes_cinereus_model.tif")
reference_layer = geotiff(SimpleSDMPredictor, joinpath("data", "input", "canada_ref_2.tif"); top=89.0)
layer_can = clip(layer, reference_layer)
layer_can = mask(replace(reference_layer, 0.0 => 1.0), layer_can)
plot(layer_can)
geotiff("data/sdms/Aeorestes_cinereus_model_can.tif", layer_can)