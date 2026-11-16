# Stationary density in the food level u. Source A has drift -delta*u and source B has
# drift delta*(1-u). At first order the local correction adds a deficit at A and a leak at B:
#
#   v_A^(1) = +(delta/2) lambda_A,  v_B^(1) = -(delta/2) lambda_B.
#
# The total current still vanishes, so rho_A = (1-u)rho and rho_B = u*rho. The total density
# is rho = f/(delta*D), f'/f = N/(delta*D), with
#
#   D = u(1-u) - (epsilon/2)[(1-u)lambda_A + u*lambda_B],
#   N = (1-u)lambda_A - u*lambda_B.

function _order(correction::Symbol)
    correction in (:leading, :local) ||
        throw(ArgumentError("correction must be :leading or :local"))
    return correction === :local ? 1.0 : 0.0
end

function _D(p::ModelParameters, c, u)
    return u * (1 - u) -
           c * p.epsilon / 2 * ((1 - u) * leave_A(p, u) + u * leave_B(p, u))
end

function _D_slope(p::ModelParameters, c, u)
    rate_slope =
        -leave_A(p, u) +
        (1 - u) * _leave_A_slope(p, u) +
        leave_B(p, u) +
        u * _leave_B_slope(p, u)
    return 1 - 2u - c * p.epsilon * rate_slope / 2
end

_N(p::ModelParameters, u) = (1 - u) * leave_A(p, u) - u * leave_B(p, u)

"""Stationary density of the leading (`:leading`) or O(epsilon) local (`:local`) PDMP."""
struct StationarySolution
    parameters::ModelParameters
    correction::Symbol
    lower::Float64
    upper::Float64
    lower_exponent::Float64
    upper_exponent::Float64
    log_normalization::Float64
    nodes::Vector{Float64}
    masses::Vector{Float64}
end

function stationary_solution(p::ModelParameters; correction::Symbol=:leading, nodes::Integer=600)
    c = _order(correction)
    D(u) = _D(p, c, u)
    lower, upper = if c * p.epsilon == 0
        (0.0, 1.0)
    else
        D(0.5) > 0 || throw(ArgumentError("epsilon too large: corrected support is empty"))
        (find_zero(D, (0.0, 0.5), Bisection()), find_zero(D, (0.5, 1.0), Bisection()))
    end
    a = _N(p, lower) / (p.delta * _D_slope(p, c, lower))
    b = _N(p, upper) / (p.delta * _D_slope(p, c, upper))
    a > 0 && b > 0 || throw(ArgumentError("non-integrable stationary density"))

    x, w = gaussjacobi(nodes, b - 1, a - 1)
    width = upper - lower
    u = @. lower + width * (1 + x) / 2
    partial = StationarySolution(p, correction, lower, upper, a, b, 0.0, u, w)
    smooth = [exp(_log_smooth(partial, ui)) for ui in u]
    masses = w .* (width / 2)^(a + b - 1) .* smooth
    total = sum(masses)
    return StationarySolution(p, correction, lower, upper, a, b, log(total), u, masses ./ total)
end

function _log_smooth(s::StationarySolution, u; points=48)
    p, c = s.parameters, _order(s.correction)
    l, h, a, b = s.lower, s.upper, s.lower_exponent, s.upper_exponent
    integrand(v) = _N(p, v) / (p.delta * _D(p, c, v)) - a / (v - l) + b / (h - v)
    middle = (l + h) / 2
    x, w = gausslegendre(points)
    log_f = (u - middle) / 2 * sum(w .* integrand.(middle .+ (u - middle) .* (x .+ 1) ./ 2))
    E = _D(p, c, u) / ((u - l) * (h - u))
    return log_f - log(p.delta * E)
end

support(s::StationarySolution) = (s.lower, s.upper)

"""Total stationary density of the food remaining at source A."""
function density(s::StationarySolution, u::Real)
    s.lower < u < s.upper || return 0.0
    l, h = s.lower, s.upper
    return exp(
        (s.lower_exponent - 1) * log(u - l) +
        (s.upper_exponent - 1) * log(h - u) +
        _log_smooth(s, u) - s.log_normalization,
    )
end

function source_density(s::StationarySolution, u::Real, source::Symbol)
    _check_source(source)
    return (source === :A ? 1 - u : u) * density(s, u)
end

"""Stationary source occupancies `(P_A, P_B)`."""
function source_probabilities(s::StationarySolution)
    return (sum(s.masses .* (1 .- s.nodes)), sum(s.masses .* s.nodes))
end

"""Stationary switching flux, equal in both directions."""
function switching_flux(s::StationarySolution)
    return sum(s.masses .* (1 .- s.nodes) .* leave_A.(Ref(s.parameters), s.nodes))
end

function cumulative(s::StationarySolution, u::Real; nodes::Integer=40)
    l, h = s.lower, s.upper
    u <= l && return 0.0
    u >= h && return 1.0
    a, b = s.lower_exponent, s.upper_exponent
    smooth(x) = exp(_log_smooth(s, x) - s.log_normalization)
    if u <= (l + h) / 2
        x, w = gaussjacobi(nodes, 0.0, a - 1)
        v = @. l + (u - l) * (1 + x) / 2
        return ((u - l) / 2)^a * sum(@. w * (h - v)^(b - 1) * smooth(v))
    end
    x, w = gaussjacobi(nodes, b - 1, 0.0)
    v = @. u + (h - u) * (1 + x) / 2
    return 1 - ((h - u) / 2)^b * sum(@. w * (v - l)^(a - 1) * smooth(v))
end

function wasserstein_distance(
    s::StationarySolution, masses::AbstractVector, food_steps::Integer
)
    length(masses) == food_steps + 1 ||
        throw(ArgumentError("masses must have food_steps + 1 entries"))
    exact = cumsum(masses)
    theory = [cumulative(s, (k + 0.5) / food_steps) for k in 0:food_steps]
    return sum(abs, theory .- exact) / food_steps
end

"""PDMP masses on the `food_steps + 1` lattice cells used by the exact chain."""
function lattice_profile(s::StationarySolution, food_steps::Integer)
    food_steps > 0 || throw(ArgumentError("food_steps must be positive"))
    edges = [cumulative(s, (k - 0.5) / food_steps) for k in 0:(food_steps + 1)]
    return diff(edges)
end
