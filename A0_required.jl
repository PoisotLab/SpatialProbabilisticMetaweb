import Pkg; Pkg.activate(".")

import CSV
using CairoMakie
using Chain
using DataFrames
using DataFramesMeta
using Distributions
using EcologicalNetworks
using EvoTrees
using GDAL
using ProgressMeter
using Random
using Shapefile
using SpeciesDistributionToolkit
using Statistics
using StatsBase
using Base.Threads: @threads
using JLD2, CodecZlib
import ColorBlendModes

include("A1_LCBD.jl")
# include("A2_colors.jl")
# include("A3_bivariate_values.jl")
# include("A4_worldshape.jl")