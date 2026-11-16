"""
    ModelParameters(; foragers=200, epsilon=0.1, delta=0.22,
                      omega0=0.02, omega1=0.2)

Parameters of the two-source ant-foraging model. `foragers` is the colony size N and
`epsilon = 1 / beta` is the inverse neutral-recruitment rate. The food level u at source A
obeys u̇ = delta(1 - u - n/N), where n foragers are at A. A lone forager leaves A or B at

    omega_A(u) = omega0 + omega1(1 - u)^2,
    omega_B(u) = omega0 + omega1 u^2.

Thus source A becomes less attractive as its food is consumed, while B becomes less
attractive as A recovers.
"""
struct ModelParameters
    foragers::Int
    epsilon::Float64
    delta::Float64
    omega0::Float64
    omega1::Float64

    function ModelParameters(
        foragers::Integer, epsilon::Real, delta::Real, omega0::Real, omega1::Real
    )
        foragers >= 2 || throw(ArgumentError("foragers must be at least two"))
        epsilon >= 0 || throw(ArgumentError("epsilon must be non-negative"))
        delta > 0 || throw(ArgumentError("delta must be positive"))
        omega0 > 0 || throw(ArgumentError("omega0 must be positive"))
        omega1 >= 0 || throw(ArgumentError("omega1 must be non-negative"))
        return new(Int(foragers), float(epsilon), float(delta), float(omega0), float(omega1))
    end
end

function ModelParameters(;
    foragers::Integer=200,
    epsilon::Real=0.1,
    delta::Real=0.22,
    omega0::Real=0.02,
    omega1::Real=0.2,
)
    return ModelParameters(foragers, epsilon, delta, omega0, omega1)
end

"""Rate at which a lone forager leaves source A, manuscript omega_+(u)."""
leave_A(p::ModelParameters, u) = p.omega0 + p.omega1 * (1 - u)^2

"""Rate at which a lone forager leaves source B, manuscript omega_-(u)."""
leave_B(p::ModelParameters, u) = p.omega0 + p.omega1 * u^2

_leave_A_slope(p::ModelParameters, u) = -2p.omega1 * (1 - u)
_leave_B_slope(p::ModelParameters, u) = 2p.omega1 * u

function _check_source(source::Symbol)
    source in (:A, :B) || throw(ArgumentError("source must be :A or :B"))
end

"""Switching rate out of `source`, where `source` is `:A` or `:B`."""
function switching_rate(p::ModelParameters, u, source::Symbol)
    _check_source(source)
    return source === :A ? leave_A(p, u) : leave_B(p, u)
end
