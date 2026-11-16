# The *reset* PDMP: the smallest process-level object that carries the complete O(ε)
# correction.
#
# Expanding the exact discrete first-order generator row out of `(+,n)` in the continuum
# gives, for a reduced test function φ,
#
#     L̄⁽¹⁾φ(+,s) = εκc₋[φ₊'(s) − φ₋'(s)] + εκ(c₋/s)[φ₋(s) − φ₊(s)]
#
# while a *local* two-state PDMP — any choice of drifts `b_i` and hazards `λ_ij` — can only
# produce the `φ₊'` half. The missing `−εκc₋φ₋'(s)` is an off-diagonal derivative, and
# `LocalPDMP` drops it. A reset PDMP generates it exactly:
#
#     (Lφ)_i = b_i(s)φ_i'(s) + λ_ij(s)[φ_j(R_ij(s)) − φ_i(s)]
#     R₊₋(s) = s(1−εκ),   R₋₊(s) = s(1+εκ)
#
# since λ⁰₊₋ = c₋/s turns φ₋(s(1−εκ)) = φ₋ − εκ s φ₋' + O(ε²) into exactly −εκc₋φ₋'.
# Physically: fixation is not instantaneous, so the population moves during the excursion
# and a switch does not land at the `s` it started from.
#
# The manuscript uses the analytic O(ε) stationary solution below.
#
# The grid above is not needed: the reset closes in elementary functions. Expanding the
# nonlocal gain term, `g(s/(1±εκ))/(1±εκ) = g ∓ εκ(sg)'`, and using `sλ⁰ = c`, the O(ε)
# stationary equations can be written as a *cross-flux*
#
#     0 = -d/ds[ b₊^ε p₊ + εκc₊p₋ ] - λ^ε₊₋p₊ + λ^ε₋₊p₋       (and the mirror)
#
# so `F₊+F₋` is constant again. Setting it to zero, the drift shifts cancel identically:
#
#     b₊^ε p₊ + εκc₊p₋ + b₋^ε p₋ - εκc₋p₊ = b₊⁰p₊ + b₋⁰p₋ = 0
#
# — the *uncorrected* drifts. Substituting `p₋ = -(b₊⁰/b₋⁰)p₊` leaves one linear ODE whose
# integrand is rational with simple poles at 0 and the two roots of
# `D = b₋⁰b₊^ε - εκc₊b₊⁰`, so partial fractions give powers again. The `1/s` residue is
# exactly 1 (`N(0)/D(0) = (1-εκβ)/(1-εκβ)`), and since `N/(sD) = O(s⁻²)` at infinity the
# residues sum to zero — hence `a₁+a₂ = 1` *exactly*, and the total density is a
# generalised arcsine, as at leading order. The branch split is independent of ε.

"""
    ResetPDMP(p::ReducedParameters)

Analytic O(ε) stationary density of the reset PDMP, with no grid or quadrature. Completes the O(ε)
generator, which `LocalPDMP`
cannot (it drops an off-diagonal derivative), while staying in closed form:

    Π₀(s) ∝ (s-ϕ₁)^{a₁-1} (ϕ₂-s)^{a₂-1},        a₁ + a₂ = 1 exactly
    Π₊/Π₀ = (1-k₂s)/(2κs),   Π₋/Π₀ = (k₁s-1)/(2κs)      (ε-free)

with `ϕ₁,ϕ₂` the roots of `D(s) = (1-k₂s)(1+εκc₋-k₁s) - εκc₊(1-k₁s)` (they agree with
`LocalPDMP`'s shifted turning points to O(ε²)) and
`a₁ = -N(ϕ₁)/(ϕ₁D'(ϕ₁))`, `N(s) = c₋(1+εκ)(1-k₂s) + c₊(1-εκ)(1-k₁s)`.

Because `a₁+a₂=1` the normalisation is `π/sin(πa₁)` and the branch weight is elementary:
`P(+) = [ϕ₁^{a₁-1}ϕ₂^{-a₁} - k₂]/(2κ)`.

# Fields
- `parameters`: the `ReducedParameters`
- `lower`, `upper`: roots `ϕ₁, ϕ₂` of `D`, the support
- `a1`, `a2`: arcsine exponents, `a1 + a2 = 1`
- `normalization`: `sin(πa₁)/π`, since `Δ^{a₁+a₂-1}B(a₁,a₂) = π/sin(πa₁)`
"""
struct ResetPDMP <: StationaryPDMP
    parameters::ReducedParameters
    lower::Float64
    upper::Float64
    a1::Float64
    a2::Float64
    normalization::Float64
end

function ResetPDMP(p::ReducedParameters)
    (; kappa, beta, epsilon) = p
    c_plus, c_minus = (1 + beta) / 2, (1 - beta) / 2
    k1, k2 = 1 + kappa, 1 - kappa
    ek = epsilon * kappa
    N(s) = c_minus * (1 + ek) * (1 - k2 * s) + c_plus * (1 - ek) * (1 - k1 * s)
    dD(s) = -k2 * (1 + ek * c_minus - k1 * s) - k1 * (1 - k2 * s) + ek * c_plus * k1
    # D(s) = A2 s² + A1 s + A0
    A2 = k1 * k2
    A1 = -(k2 * (1 + ek * c_minus) + k1) + ek * c_plus * k1
    A0 = (1 + ek * c_minus) - ek * c_plus
    disc = sqrt(A1^2 - 4 * A2 * A0)
    lower, upper = minmax((-A1 - disc) / (2A2), (-A1 + disc) / (2A2))
    a1 = -N(lower) / (lower * dD(lower))
    a2 = -N(upper) / (upper * dD(upper))
    return ResetPDMP(p, lower, upper, a1, a2, sin(π * a1) / π)
end

support(s::ResetPDMP) = (s.lower, s.upper)

_total(s::ResetPDMP, x) = s.normalization * (x - s.lower)^(s.a1 - 1) * (s.upper - x)^(s.a2 - 1)

function branch_density(s::ResetPDMP, x::Real, branch::Symbol)
    _check_branch(branch)
    (; kappa) = s.parameters
    k1, k2 = 1 + kappa, 1 - kappa
    split = branch === :plus ? (1 - k2 * x) : (k1 * x - 1)
    return _total(s, x) * split / (2 * kappa * x)
end

# Π₀ is exactly Beta(a₁, a₂) on (ϕ₁, ϕ₂) with a₁+a₂ = 1, and the branch split is
# `Π₊/Π₀ = (1-k₂s)/(2κs)`, `Π₋/Π₀ = (k₁s-1)/(2κs)`. So everything needs only
#
#     E[1/s]  = ϕ₁^{a₁-1} ϕ₂^{-a₁},
#     E[1/s²] = E[1/s] (a₂/ϕ₁ + a₁/ϕ₂),
#
# the first from `₂F₁(1,a₁;1;z) = (1-z)^{-a₁}`, the second by differentiating the first with
# respect to a common shift of ϕ₁ and ϕ₂.

function _inverse_moments(s::ResetPDMP)
    inv1 = s.lower^(s.a1 - 1) * s.upper^(-s.a1)
    return inv1, inv1 * (s.a2 / s.lower + s.a1 / s.upper)
end

function branch_probabilities(s::ResetPDMP)
    (; kappa) = s.parameters
    inv1, _ = _inverse_moments(s)
    return ((inv1 - (1 - kappa)) / (2 * kappa), ((1 + kappa) - inv1) / (2 * kappa))
end

"""
    switching_flux(s::ResetPDMP) -> Float64

Mean of `∫ λ₊Π₊ ds` and `∫ λ₋Π₋ ds`, with the corrected hazards. The analytic reset solution
is complete only to O(ε), so the two directions agree to O(ε²) rather than exactly.
"""
function switching_flux(s::ResetPDMP)
    (; kappa, beta, epsilon) = s.parameters
    k1, k2 = 1 + kappa, 1 - kappa
    c_plus, c_minus = (1 + beta) / 2, (1 - beta) / 2
    ek = epsilon * kappa
    inv1, inv2 = _inverse_moments(s)
    plus = c_minus * (1 + ek) * (inv2 - k2 * inv1) / (2 * kappa)
    minus = c_plus * (1 - ek) * (k1 * inv1 - inv2) / (2 * kappa)
    return (plus + minus) / 2
end

stationary_moments(s::ResetPDMP) = _beta_moments(s.lower, s.upper, s.a1, s.a2)
