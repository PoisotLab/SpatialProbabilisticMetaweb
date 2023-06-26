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
using GLMakie
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

include("A1_LCBD.jl")
include("A2_colors.jl")
include("A3_bivariate_values.jl")
include("A4_worldshape.jl")