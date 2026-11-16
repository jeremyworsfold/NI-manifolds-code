module Ants

# Foraging ants choosing between sources A and B. Neutral recruitment is fast, while the
# food remaining at source A changes slowly and controls the rates at which lone foragers
# change source. The reduced process is a two-source PDMP in the food level u.
using LinearAlgebra
using Random
using SparseArrays
using FastGaussQuadrature: gaussjacobi, gausslegendre
using Roots: find_zero, Bisection
using OrdinaryDiffEqTsit5: ODEProblem, Tsit5, solve, DiscreteCallback, terminate!
using OrdinaryDiffEqSDIRK: TRBDF2, MatrixOperator
using LinearSolve: KLUFactorization
using JumpProcesses: JumpProblem, VariableRateJump, VR_Direct
using Statistics: mean
using StatsBase: ecdf
using ..NIManifoldsCode: stationary_vector

export ModelParameters,
    leave_A,
    leave_B,
    switching_rate,
    StationarySolution,
    stationary_solution,
    support,
    density,
    source_density,
    source_probabilities,
    switching_flux,
    cumulative,
    wasserstein_distance,
    lattice_profile,
    survival,
    residence_moments,
    transit_time,
    merged_residences,
    sample_moments,
    empirical_survival,
    ExactChain,
    exact_chain,
    food_masses,
    exact_residence

include("model.jl")
include("theory/stationary.jl")
include("theory/residence.jl")
include("theory/merged.jl")
include("exact.jl")

end
