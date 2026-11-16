# Stationary laws of the reduced TK2 PDMP, in the slow coordinate s = n/V of the total
# population. There is one type per correction of the manuscript:
#
#     LeadingPDMP   stationary/leading.jl   instantaneous fixation, a Beta law
#     LocalPDMP     stationary/local.jl     first-order excursion correction, local part only
#     ResetPDMP     stationary/reset.jl     complete first-order correction, with the reset
#
# Each file holds the type's density, normalisation and observables (branch probabilities,
# switching flux, moments). This file holds what they share: the abstract type, the entry point,
# and the two integrals that reduce those observables to elementary functions.

"""
Stationary law of a reduced TK2 PDMP on the branches `:plus` and `:minus`. Concrete types are
`LeadingPDMP`, `LocalPDMP` and `ResetPDMP`. All of them store their `ReducedParameters` as
`parameters` and implement `support`, `branch_density`, `branch_probabilities`,
`switching_flux` and `stationary_moments`.
"""
abstract type StationaryPDMP end

"""
    stationary_solution(p::ReducedParameters; correction = :leading)

Stationary law of the leading (`:leading`), local (`:local`) or reset (`:reset`) PDMP.

The corrections are proportional to κ, so at κ = 0 every order is the leading one. The
corrected formulas divide by κ and lose all precision as κ → 0, so for |κ| < 1e-6, where the
O(εκ) correction is negligible, the leading law is returned for every `correction`.
"""
function stationary_solution(p::ReducedParameters; correction::Symbol=:leading)
    correction in (:leading, :local, :reset) ||
        throw(ArgumentError("correction must be :leading, :local or :reset"))
    (correction === :leading || abs(p.kappa) < 1e-6) && return LeadingPDMP(p)
    return correction === :local ? LocalPDMP(p) : ResetPDMP(p)
end

"""Total stationary density `Π₊(s) + Π₋(s)` of the slow coordinate."""
density(s::StationaryPDMP, x::Real) = branch_density(s, x, :plus) + branch_density(s, x, :minus)

function _check_branch(branch::Symbol)
    branch in (:plus, :minus) || throw(ArgumentError("branch must be :plus or :minus"))
end

# --- Euler integrals -----------------------------------------------------------------
#
# All three stationary densities are products of powers on the support (ϕ₁, ϕ₂),
#
#     (s-ϕ₁)^{α-1} (ϕ₂-s)^{β-1} s^{-γ},
#
# so every observable is an Euler integral. With s = ϕ₁ + Δu, Δ = ϕ₂ - ϕ₁, and a Pfaff
# transformation to put the argument in (0,1),
#
#     ∫ (s-ϕ₁)^{α-1}(ϕ₂-s)^{β-1} s^{-γ} ds
#         = Δ^{α+β-1} B(α,β) ϕ₁^{α-γ} ϕ₂^{-α} ₂F₁(α+β-γ, α; α+β; Δ/ϕ₂).
#
# The leading and reset densities are Beta laws (α+β = 1, γ = 0 for the total), so their branch
# weights and moments collapse further to elementary closed forms. Nothing samples the density
# on a grid, so the endpoint singularities cost no accuracy.

"""`∫_{ϕ₁}^{ϕ₂} (s-ϕ₁)^{α-1}(ϕ₂-s)^{β-1} s^{-γ} ds`, in closed form (see the comment above)."""
function _euler_integral(ϕ1, ϕ2, α, β, γ)
    Δ = ϕ2 - ϕ1
    return Δ^(α + β - 1) *
           SpecialFunctions.beta(α, β) *
           ϕ1^(α - γ) *
           ϕ2^(-α) *
           _₂F₁(α + β - γ, α, α + β, Δ / ϕ2)
end

"""Mean, second moment and variance of a Beta(a, b) law pushed onto `(ϕ1, ϕ2)`."""
function _beta_moments(ϕ1, ϕ2, a, b)
    Δ = ϕ2 - ϕ1
    mean = ϕ1 + Δ * a / (a + b)
    variance = Δ^2 * a * b / ((a + b)^2 * (a + b + 1))
    return (; mean, second=variance + mean^2, variance)
end

"""
Mean branch coordinate `y = P₊ - P₋` of the stationary PDMP, the counterpart of
`mean_branch_coordinate(::ExactChain)`.
"""
function mean_branch_coordinate(s::StationaryPDMP)
    plus, minus = branch_probabilities(s)
    return plus - minus
end
