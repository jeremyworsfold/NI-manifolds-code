"""
    LeadingPDMP(p::ReducedParameters)

Leading-order stationary law of the PDMP, where fixation is instantaneous. The total density
of `s` is a Beta(c₋, c₊) law on `(1/k₁, 1/k₂)`, with `c± = (1 ± β)/2` and branch decay rates
`k₁ = 1 + κ`, `k₂ = 1 - κ`, split between the two branches. It does not depend on ε.

# Fields
- `parameters`: the `ReducedParameters` (only κ and β enter)
- `k1`, `k2`: branch decay rates `1 ± κ`
- `normalization`: constant `C` in `Π±(s) = C (k₁s-1)^{…}(1-k₂s)^{…}/s`
"""
struct LeadingPDMP <: StationaryPDMP
    parameters::ReducedParameters
    k1::Float64
    k2::Float64
    normalization::Float64
end

function LeadingPDMP(p::ReducedParameters)
    (; kappa, beta) = p
    normalization =
        cos(π * beta / 2) * (1 + kappa) * ((1 - kappa) / (1 + kappa))^((1 - beta) / 2) /
        (2π * kappa)
    return LeadingPDMP(p, 1 + kappa, 1 - kappa, normalization)
end

"""Support of the stationary density, the two deterministic branch fixed points `(1/k₁, 1/k₂)`."""
support(s::LeadingPDMP) = (1 / s.k1, 1 / s.k2)

function branch_density(s::LeadingPDMP, x::Real, branch::Symbol)
    _check_branch(branch)
    (; k1, k2, normalization) = s
    expo = branch === :plus ? -(1 + s.parameters.beta) / 2 : (1 - s.parameters.beta) / 2
    return (k1 * x - 1)^expo * (1 - k2 * x)^(-expo) * normalization / x
end

"""Moments of `s` under the stationary density, which is exactly Beta(c₋, c₊) on the support."""
function stationary_moments(s::LeadingPDMP)
    beta = s.parameters.beta
    c_plus, c_minus = (1 + beta) / 2, (1 - beta) / 2
    return _beta_moments(support(s)..., c_minus, c_plus)
end

"""Stationary probabilities `(P₊, P₋)` of occupying the plus and minus branches."""
function branch_probabilities(s::LeadingPDMP)
    (; kappa, beta) = s.parameters
    c_plus, c_minus = (1 + beta) / 2, (1 - beta) / 2
    if abs(kappa) < 1e-6
        # The removable κ = 0 limit, with the first-order terms retained for continuity.
        return (c_plus * (1 - c_minus * kappa), c_minus * (1 + c_plus * kappa))
    end
    r = (1 + kappa) / (1 - kappa)
    P_plus = (1 - kappa) / (2kappa) * (r^c_plus - 1)
    P_minus = (1 + kappa) / (2kappa) * (1 - r^(-c_minus))
    return (P_plus, P_minus)
end

"""
    switching_flux(s::LeadingPDMP) -> Float64

Stationary number of switches per unit time in either direction,
`c₊c₋(1-κ)^c₋(1+κ)^c₊`.
"""
function switching_flux(s::LeadingPDMP)
    (; kappa, beta) = s.parameters
    c_plus, c_minus = (1 + beta) / 2, (1 - beta) / 2
    return c_plus * c_minus * (1 - kappa)^c_minus * (1 + kappa)^c_plus
end

"""
    sample_residence_time(s::LeadingPDMP, x, branch, u) -> (t, x_exit)

Sample one leading-order residence exactly by inverting its closed-form survival law, as in the
schematic path of the manuscript. The PDMP enters `branch` (`:plus` or `:minus`) at `x`; `u`
lies in `(0, 1]`. Returns the residence time and the state at which the branch is left.
"""
function sample_residence_time(s::LeadingPDMP, x::Real, branch::Symbol, u::Real)
    _check_branch(branch)
    beta = s.parameters.beta
    c_plus, c_minus = (1 + beta) / 2, (1 - beta) / 2
    k, c = branch === :plus ? (s.k1, c_minus) : (s.k2, c_plus)
    A = k * x * (u^(-1 / c) - 1)
    fixed_point = 1 / k
    return (log1p(A) / k, fixed_point + (x - fixed_point) / (1 + A))
end

"""
One residence of a leading-order PDMP path: times `t` and states `x` along the deterministic
flow on `branch` (`:plus` or `:minus`).
"""
struct PathSegment
    t::Vector{Float64}
    x::Vector{Float64}
    branch::Symbol
end

"""
    sample_path(s::LeadingPDMP, x0, residences; rng, points = 60) -> Vector{PathSegment}

Exact sample path of the leading-order PDMP, entering the plus branch at `x0`, with
`residences` residences that alternate between the branches. Each residence is drawn by
`sample_residence_time`, so there is no time-stepping error. Between switches the state relaxes
to the branch fixed point `1/k` as `x(t) = 1/k + (x - 1/k) e^{-kt}`, sampled at `points` times.
"""
function sample_path(s::LeadingPDMP, x0::Real, residences::Integer; rng, points::Integer=60)
    segments = PathSegment[]
    x, t, branch = Float64(x0), 0.0, :plus
    for _ in 1:residences
        tau, x_exit = sample_residence_time(s, x, branch, 1 - rand(rng))
        k = branch === :plus ? s.k1 : s.k2
        times = range(0, tau; length=points)
        push!(
            segments,
            PathSegment(collect(t .+ times), @.(1 / k + (x - 1 / k) * exp(-k * times)), branch),
        )
        x, t = x_exit, t + tau
        branch = branch === :plus ? :minus : :plus
    end
    return segments
end
