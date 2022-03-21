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