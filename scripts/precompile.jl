# Copyright 2023, Ben Cardoen
# This script is used to generate a precompilation trace
# It'll invoke the tests, so any code tested is precompiled subsequently.
using ArgParse
using SmlmTools, Images
using CSV, DataFrames
using Logging, LoggingExtras, Dates
include(joinpath(pkgdir(SmlmTools), "test", "runtests.jl"))
