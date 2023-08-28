function background_map(
    fig = Figure();
    lims=(left = -142.0, right = -52.0, bottom = 41.0, top = 84.0),
    tuple=false,
    kw...
)
    shapes = Shapefile.shapes(Shapefile.Table("data/shapefiles/land/land_50m_curved.shp"))
    ga = GeoAxis(
        fig isa Figure ? fig[1,1] : fig;
        source = "+proj=longlat +datum=WGS84",
        dest = "esri:102002", # Lambert Conformal Conic
        lonlims = (lims.left, lims.right),
        latlims = (lims.bottom, lims.top),
        xlabel = "Longitude",
        ylabel = "Latitude",
        kw...
    )
    foreach(shapes) do sh
        poly!(ga, sh; shading=false, strokecolor=:darkgrey, strokewidth=1, color=:transparent)
    end
    if tuple
        return (fig=fig, ga=ga)
    else
        return fig
    end
end

function Makie.convert_arguments(::Type{<:Poly}, p::Shapefile.Polygon)
    polys = Shapefile.GeoInterface.coordinates(p)
    ps = map(polys) do pol
        Polygon(
            Point2f.(pol[1]), # interior
            map(x -> Point2f.(x), pol[2:end])
        )
    end
    (ps,)
end