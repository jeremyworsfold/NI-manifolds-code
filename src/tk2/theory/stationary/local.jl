"""
    LocalPDMP(p::ReducedParameters)

The first-order excursion-interruption correction to the symmetric PDMP, in the
continuum limit `V→∞` at fixed `ε = D·V`. It has two pieces, *both* proportional to the branch-decay
asymmetry κ and both vanishing identically at `k₊ = k₋`:

- a constant drift shift `b₊(s) = 1 + ε·κ·c₋ - k₁·s`, `b₋(s) = 1 - ε·κ·c₊ - k₂·s` — the
  population decays at a *mixture* of `k₊` and `k₋` while an excursion is in flight, so
  each branch is pulled toward the interior. This moves the deterministic turning points,
  hence the `support` method below, and it is the dominant of the two;
- a switching-hazard shift `λ(s) = (c/s)(1 ± ε·κ)`.

The stationary density is **closed form**, like `LeadingPDMP`'s: the correction
moves the poles and residues of the J=0
integrand but leaves it rational, so no quadrature is needed (see the "Closed form"
section below). It plugs into
`finite_size_profile` unchanged. At `ε = 0` it reduces exactly to `LeadingPDMP`.

Accuracy: at the reduced-chain level the correction is first-order complete — the residual
against the exact chain's symmetric fast-fixation projection is `O(ε²)`, with the corrected/leading error
ratio numerically equal to `ε`. Most of *this*
object's own error is not the reduction but the continuum+LNA layer built on top of it: at `V=30` it
is essentially all of it, at `V=100, ε=0.2` about half. The finite-volume reduced chain
is the more accurate object at modest `V`; this one is the closed-form, interpretable one.

# Fields
- `parameters`: the `ReducedParameters` (κ, β and the controlling parameter ε = D·V)
- `lower`, `upper`: the shifted turning points `ϕ₁, ϕ₂`, the support
- `a1`, `a2`: residues of the J=0 integrand at `ϕ₁` and `ϕ₂`
- `normalization`: constant `C` of the closed-form density
"""
struct LocalPDMP <: StationaryPDMP
    parameters::ReducedParameters
    lower::Float64
    upper::Float64
    a1::Float64
    a2::Float64
    normalization::Float64
end

# --- Closed form -------------------------------------------------------------------
#
# The correction preserves both structural features that make the leading-order solution
# elementary: the hazards keep their `c/s` form (only the prefactors move,
# `c∓ → c∓(1±εκ)`) and the drifts stay linear in `s` (only the intercept moves). So the
# J=0 integrand is still rational with simple poles at 0, ϕ₁, ϕ₂, and partial fractions give
#
#   d(log u)/ds = a₁/(s-ϕ₁) - (a₁+a₂)/s - a₂/(ϕ₂-s)
#
# with the *residues* shifted from `c∓` to `a₁,a₂` and the poles from `1/k` to `ϕ`. Hence
# `Π₊ ∝ (s-ϕ₁)^{a₁-1}(ϕ₂-s)^{a₂} s^{-P}/k₁` and `Π₋ ∝ (s-ϕ₁)^{a₁}(ϕ₂-s)^{a₂-1} s^{-P}/k₂`,
# `P = a₁+a₂`.

function LocalPDMP(p::ReducedParameters)
    (; kappa, beta, epsilon) = p
    c_plus, c_minus = (1 + beta) / 2, (1 - beta) / 2
    k1, k2 = 1 + kappa, 1 - kappa
    ek = epsilon * kappa
    # The deterministic turning points, shifted by the first-order excursion drift: the roots
    # of `b±(s) = 0` with the constant shift included. Both move inward (toward each other)
    # for κ > 0.
    lower, upper = (1 + ek * c_minus) / k1, (1 - ek * c_plus) / k2
    # Residues of the J=0 integrand at `ϕ₁` and `ϕ₂`: `a₁ = c₋(1+εκ)/(1+εκ·c₋)`,
    # `a₂ = c₊(1-εκ)/(1-εκ·c₊)`, reducing to `c₋, c₊` at ε=0. Their sum is
    # `P = 1 - ε²κ²c₊c₋ + O(ε³)` — the departure from 1 is second order, which is what makes
    # the normalisation below collapse.
    a1 = c_minus * (1 + ek) / (1 + ek * c_minus)
    a2 = c_plus * (1 - ek) / (1 - ek * c_plus)
    # Normalisation of the closed form. Summing the two branches,
    # `Π₀ ∝ κ(2s-ε) s^{-P} (s-ϕ₁)^{a₁-1}(ϕ₂-s)^{a₂-1}/(k₁k₂)`, and both Euler integrals
    # reduce: the `ε` term degenerates exactly (its `₂F₁` has first parameter equal to
    # third), and the `2s` term leaves `₂F₁(P-1,a₁;P;-Δ/ϕ₁)`, which a Pfaff transformation
    # turns into `(ϕ₁/ϕ₂)^{a₁}·₂F₁(1,a₁;P;Δ/ϕ₂)`, with argument `Δ/ϕ₂ = 1-ϕ₁/ϕ₂ ∈ (0,1)`. So
    #
    #     1/C = κ/(k₁k₂) · Δ^{P-1} · B(a₁,a₂) · ϕ₁^{-a₂}ϕ₂^{-a₁} · (2ϕ₁·₂F₁(1,a₁;P;Δ/ϕ₂) - ε).
    #
    # At ε=0 the `₂F₁` collapses to `(ϕ₂/ϕ₁)^{c₋}` and this is `2κ·B(c₋,c₊)/(k₁k₂)`, i.e.
    # the leading-order constant. Agrees with adaptive quadrature of the density to ~1e-11.
    P = a1 + a2
    Δ = upper - lower
    F = _₂F₁(1, a1, P, Δ / upper)
    inverse_normalization =
        kappa / (k1 * k2) *
        Δ^(P - 1) *
        SpecialFunctions.beta(a1, a2) *
        lower^(-a2) *
        upper^(-a1) *
        (2 * lower * F - epsilon)
    return LocalPDMP(p, lower, upper, a1, a2, 1 / inverse_normalization)
end

support(s::LocalPDMP) = (s.lower, s.upper)

"""Unnormalised branch densities `(Π₊, Π₋)` at `x`, from the closed form above."""
function _shape(s::LocalPDMP, x)
    (; lower, upper, a1, a2) = s
    k1, k2 = 1 + s.parameters.kappa, 1 - s.parameters.kappa
    base = x^(-(a1 + a2))
    return (x - lower)^(a1 - 1) * (upper - x)^a2 * base / k1,
    (x - lower)^a1 * (upper - x)^(a2 - 1) * base / k2
end

function branch_density(s::LocalPDMP, x::Real, branch::Symbol)
    _check_branch(branch)
    return s.normalization * _shape(s, x)[branch === :plus ? 1 : 2]
end

# Multiplying `Π_branch` by `s^m` shifts `γ` to `P-m`, and the hazards `c∓(1±εκ)/s` shift it
# to `P+1`, where the ₂F₁ has first parameter zero and equals 1.

"""`∫ s^m Π_branch ds` for the locally corrected density."""
function _branch_integral(s::LocalPDMP, branch::Symbol, m)
    (; lower, upper, a1, a2) = s
    P = a1 + a2
    α, β, k = if branch === :plus
        (a1, a2 + 1, 1 + s.parameters.kappa)
    else
        (a1 + 1, a2, 1 - s.parameters.kappa)
    end
    return s.normalization * _euler_integral(lower, upper, α, β, P - m) / k
end

function branch_probabilities(s::LocalPDMP)
    plus, minus = _branch_integral(s, :plus, 0), _branch_integral(s, :minus, 0)
    return (plus / (plus + minus), minus / (plus + minus))
end

"""
    switching_flux(s::LocalPDMP) -> Float64

`J = ∫ λ₊Π₊ ds = ∫ λ₋Π₋ ds` with the corrected hazards `λ₊ = c₋(1+εκ)/s`,
`λ₋ = c₊(1-εκ)/s`. The two are equal because the zero-current solution is exact for the
local PDMP; their mean is returned.
"""
function switching_flux(s::LocalPDMP)
    (; kappa, beta, epsilon) = s.parameters
    c_plus, c_minus = (1 + beta) / 2, (1 - beta) / 2
    ek = epsilon * kappa
    plus = c_minus * (1 + ek) * _branch_integral(s, :plus, -1)
    minus = c_plus * (1 - ek) * _branch_integral(s, :minus, -1)
    return (plus + minus) / 2
end

function stationary_moments(s::LocalPDMP)
    total(m) = _branch_integral(s, :plus, m) + _branch_integral(s, :minus, m)
    norm = total(0)
    mean = total(1) / norm
    second = total(2) / norm
    return (; mean, second, variance=second - mean^2)
end
