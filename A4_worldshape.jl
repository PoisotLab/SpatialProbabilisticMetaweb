function _download_shapefile(res)
    @assert res ∈ [50, 100]
    fn = "ne_$(res)m_land.shp"
    ispath("shapefiles") || mkpath("shapefiles")
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

function plot(layer::T, ws::Shapefile.Handle; ws_c=:lightgrey, ws_lc=:grey, ws_lw=0.5, kw...) where {T<:SimpleSDMLayer}
    plot(ws; c=ws_c, lc=ws_lc, lw=ws_lw)
    plot!(layer; kw...)
end

function bivariate(l1::T1, l2::T2, ws::Shapefile.Handle; ws_c=:lightgrey, ws_lc=:grey, ws_lw=0.5, kw...) where {T1 <: SimpleSDMLayer, T2 <: SimpleSDMLayer}
    plot(ws; c=ws_c, lc=ws_lc, lw=ws_lw, cb=:none)
    bivariate!(l1, l2; kw...)
end