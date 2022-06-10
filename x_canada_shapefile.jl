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

coords = (left=-145.0, right=-50.0, bottom=40.0)
layer = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
@time layer_mask = mask(canada_poly, layer) # 140 sec.
plot(layer_mask)