import Pkg; Pkg.activate(".")

using ArchGDAL
import CSV
using CairoMakie
import ColorBlendModes
using Chain
using DataFrames
using DataFramesMeta
using Distributions
using Downloads
using EcologicalNetworks
using EvoTrees
using GDAL
using GeoMakie
if get(ENV, "CC_CLUSTER", "local") == "local"
    # Load GLMakie only if local, not on clusters
    using CairoMakie
end
using JLD2, CodecZlib
using ProgressMeter
using Random
using Shapefile
using SpeciesDistributionToolkit
using Statistics
using StatsBase
using Base.Threads: @threads
using ZipFile

using SpeciesDistributionToolkit: _read_geotiff as read_geotiff
using SpeciesDistributionToolkit: _write_geotiff as write_geotiff
using SpeciesDistributionToolkit: boundingbox

include("scripts/lib/A1_LCBD.jl")
include("scripts/lib/A2_colors.jl")
include("scripts/lib/A3_bivariate_values.jl")
include("scripts/lib/A4_worldshape.jl")

# Set default resolution for figures
set_theme!()
update_theme!(
    CairoMakie=(; px_per_unit=3),
)
