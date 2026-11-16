module TK2

using Catalyst
using FiniteStateProjection: FSPSystem, SteadyState
using Integrals
using SpecialFunctions
using SparseArrays
using HypergeometricFunctions: _₂F₁
using ..NIManifoldsCode: stationary_vector

# `reduce` is deliberately not exported: it would clash with `Base.reduce` for anyone who
# writes `using NIManifoldsCode.TK2`. Call it as `TK2.reduce`.
export ModelParameters,
    ReducedParameters,
    realize,
    StationaryPDMP,
    LeadingPDMP,
    LocalPDMP,
    ResetPDMP,
    stationary_solution,
    support,
    density,
    branch_density,
    branch_probabilities,
    switching_flux,
    stationary_moments,
    sample_residence_time,
    PathSegment,
    sample_path,
    finite_size_profile,
    density_l1,
    density_w1,
    ExactChain,
    exact_chain,
    default_cutoff,
    population_masses,
    boundary_mass,
    mean_branch_coordinate

include("model.jl")
include("theory/stationary.jl")
include("theory/stationary/leading.jl")
include("theory/stationary/local.jl")
include("theory/stationary/reset.jl")
include("theory/finite_size.jl")
include("exact.jl")

end
