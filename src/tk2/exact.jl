# Exact microscopic chain on the populations (n₁, n₂) of the two species, truncated to the
# square 0 ≤ n₁, n₂ ≤ cutoff by finite state projection (FiniteStateProjection drops the
# transitions that leave the square). State (n₁, n₂) has index n₁ + (cutoff + 1) n₂ + 1.

"""Truncated exact chain with its stationary law (generator rows = from-state)."""
struct ExactChain
    parameters::ModelParameters
    cutoff::Int
    generator::SparseMatrixCSC{Float64,Int}
    stationary::Vector{Float64}
end

"""
    default_cutoff(p::ModelParameters; sigma = 8, extra = 20)

Population cutoff per species: `sigma` Poisson widths beyond the upper branch fixed point
`N/(1 - |κ|)`, where the stationary law of the total population ends, plus `extra`.
"""
function default_cutoff(p::ModelParameters; sigma::Real=8, extra::Integer=20)
    upper = p.population / (1 - abs(p.decay_bias))
    return ceil(Int, upper + sigma * sqrt(upper)) + extra
end

"""Rectangular finite state projection of the microscopic TK2 reaction system, solved for its stationary law."""
function exact_chain(p::ModelParameters; cutoff::Integer=default_cutoff(p))
    cutoff >= 1 || throw(ArgumentError("cutoff must be positive"))
    system, parameters = _reaction_system(p)
    forward = SparseMatrixCSC(
        FSPSystem(system), (cutoff + 1, cutoff + 1), parameters, SteadyState()
    )
    generator = sparse(transpose(forward))
    return ExactChain(p, cutoff, generator, stationary_vector(generator))
end

"""States `(n₁, n₂)` in the order of the chain's indices."""
_states(chain::ExactChain) = [(n1, n2) for n2 in 0:chain.cutoff for n1 in 0:chain.cutoff]

"""Stationary masses of the total population n = n₁ + n₂ = 0, …, 2·cutoff."""
function population_masses(chain::ExactChain)
    masses = zeros(2chain.cutoff + 1)
    for ((n1, n2), weight) in zip(_states(chain), chain.stationary)
        masses[n1 + n2 + 1] += weight
    end
    return masses
end

"""Stationary mass of the states on the cutoff, which the truncation discards flux from."""
function boundary_mass(chain::ExactChain)
    sum(
        weight for
        ((n1, n2), weight) in zip(_states(chain), chain.stationary) if max(n1, n2) == chain.cutoff
    )
end

# The symmetric fast-fixation projection assigns state (n₁, n₂) to the plus branch with
# probability q₊ = n₁/n, n = n₁ + n₂ (the committor), and to the minus branch with 1 - q₊. The
# empty state has no branch and is split evenly. This is the observable the theory predicts:
# compare against it, never against a raw marginal.

"""Stationary branch occupancies `(P₊, P₋)` under the symmetric fast-fixation projection."""
function branch_probabilities(chain::ExactChain)
    plus = sum(
        weight * (n1 + n2 == 0 ? 0.5 : n1 / (n1 + n2)) for
        ((n1, n2), weight) in zip(_states(chain), chain.stationary)
    )
    return (plus, 1 - plus)
end

"""
Mean branch coordinate `y = E[(n₁ - n₂)/(n₁ + n₂)]` over the occupied states `n₁ + n₂ > 0`,
which is `P₊ - P₋` of the symmetric fast-fixation projection restricted to them.
"""
function mean_branch_coordinate(chain::ExactChain)
    occupied, total = 0.0, 0.0
    for ((n1, n2), weight) in zip(_states(chain), chain.stationary)
        n1 + n2 == 0 && continue
        occupied += weight
        total += weight * (n1 - n2) / (n1 + n2)
    end
    return total / occupied
end
