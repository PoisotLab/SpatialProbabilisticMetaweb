using Downloads
using ZipFile: ZipFile

ispath("ecoregions") || mkpath("ecoregions")

for scale in ["district", "region", "province", "zone"]
    _url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/$(scale)/eco$(scale)_shp.zip"
    _tmp = Downloads.download(_url)
    _zip = ZipFile.Reader(_tmp)
    for _zf in _zip.files
        println("Filename: $(_zf.name)")
        tname = last(split(_zf.name, "/"))
        if !isempty(tname)
            open(joinpath("ecoregions", tname), "w") do io
                write(io, read(_zf, String))
            end
        end
    end
    close(_zip)
end
