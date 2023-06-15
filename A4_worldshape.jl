function _download_shapefile(res)
    @assert res ∈ [50, 100]
    fn = "ne_$(res)m_land.shp"
    ispath(joinpath("shapefiles", "land")) || mkpath(joinpath("shapefiles", "land"))
    tf = joinpath("shapefiles", fn)
    if !isfile(tf)
        dir = "https://github.com/nvkelso/natural-earth-vector/" *
            "raw/master/$(res)m_physical/"
        download("$dir/$fn", tf)
    end
    return tf
end

function worldshape(res)
    handle = open(_download_shapefile(res), "r") do io
        read(io, Shapefile.Handle)
    end
    return handle
end

# function plot(layer::T, ws::Shapefile.Handle; ws_c=:lightgrey, ws_lc=:grey, ws_lw=0.5, kw...) where {T<:SimpleSDMLayer}
#     plot(ws; c=ws_c, lc=ws_lc, lw=ws_lw, xaxis="Longitude", yaxis="Latitude")
#     plot!(layer; kw...)
# end

# function bivariate(l1::T1, l2::T2, ws::Shapefile.Handle; ws_c=:lightgrey, ws_lc=:grey, ws_lw=0.5, kw...) where {T1 <: SimpleSDMLayer, T2 <: SimpleSDMLayer}
#     plot(ws; c=ws_c, lc=ws_lc, lw=ws_lw, cb=:none, xaxis="Longitude", yaxis="Latitude")
#     bivariate!(l1, l2; kw...)
# end

function background_map(;
    fig=nothing,
    pos=[1,1],
    lims=(left=-145.0, right=-50.0, bottom=40.0, top=89.0),
    kw...
)
    shapes = Shapefile.shapes(Shapefile.Table("shapefiles/land/land_50m_curved.shp"))
    if fig==nothing
        fig = Figure()
    end
    ga = GeoAxis(
        fig[pos...];
        source = "+proj=longlat +datum=WGS84",
        dest = "esri:102002", # Lambert Conformal Conic
        lonlims = (lims.left, lims.right),
        latlims = (lims.bottom, lims.top),
        xlabel = "Longitude",
        ylabel = "Latitude",
        kw...
    )
    foreach(shapes) do sh
        poly!(ga, sh; shading=false, strokecolor=:darkgrey, strokewidth=1, color=:lightgrey)
    end
    return fig
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