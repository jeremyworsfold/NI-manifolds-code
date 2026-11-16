# Shared by every example script.

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(@__DIR__, "run.jl"))
include(joinpath(@__DIR__, "comparison.jl"))
include(joinpath(@__DIR__, "figure_data.jl"))
