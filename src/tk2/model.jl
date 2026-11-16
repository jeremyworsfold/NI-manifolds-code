"""Microscopic parameters for the two-species Togashi--Kaneko model."""
struct ModelParameters
    population::Int
    diffusion::Float64
    decay_bias::Float64
    influx_bias::Float64

    function ModelParameters(
        population::Integer, diffusion::Real, decay_bias::Real, influx_bias::Real
    )
        population > 0 || throw(ArgumentError("population must be positive"))
        diffusion >= 0 || throw(ArgumentError("diffusion must be non-negative"))
        abs(decay_bias) < 1 || throw(ArgumentError("abs(decay_bias) must be < 1"))
        abs(influx_bias) <= 1 || throw(ArgumentError("abs(influx_bias) must be <= 1"))
        return new(Int(population), float(diffusion), float(decay_bias), float(influx_bias))
    end
end

function ModelParameters(;
    population::Integer=100, diffusion::Real=1e-3, decay_bias::Real=0.5, influx_bias::Real=0.0
)
    ModelParameters(population, diffusion, decay_bias, influx_bias)
end

"""Dimensionless parameters used by TK2's reduced descriptions."""
struct ReducedParameters
    kappa::Float64
    beta::Float64
    epsilon::Float64

    function ReducedParameters(kappa::Real, beta::Real, epsilon::Real)
        0 <= kappa < 1 || throw(ArgumentError("need 0 <= kappa < 1"))
        abs(beta) <= 1 || throw(ArgumentError("abs(beta) must be <= 1"))
        epsilon >= 0 || throw(ArgumentError("epsilon must be non-negative"))
        return new(float(kappa), float(beta), float(epsilon))
    end
end

function ReducedParameters(; kappa::Real=0.5, beta::Real=0.0, epsilon::Real=0.1)
    ReducedParameters(kappa, beta, epsilon)
end

"""Convert microscopic TK2 parameters to the dimensionless reduced parameters."""
function reduce(p::ModelParameters)
    ReducedParameters(p.decay_bias, p.influx_bias, p.diffusion * p.population)
end

"""Realise reduced TK2 parameters at a chosen characteristic population."""
function realize(p::ReducedParameters; population::Integer=100)
    return ModelParameters(population, p.epsilon / population, p.kappa, p.beta)
end

"""Catalyst reaction system and parameter map for the microscopic model."""
function _reaction_system(p::ModelParameters)
    N = p.population
    k1 = p.diffusion * (1 + p.decay_bias)
    k2 = p.diffusion * (1 - p.decay_bias)
    b1 = p.diffusion * (1 + p.influx_bias) * N / 2
    b2 = p.diffusion * (1 - p.influx_bias) * N / 2
    r = 1 / N
    system = @reaction_network begin
        @species X₁(t) X₂(t)
        k1, X₁ --> ∅
        k2, X₂ --> ∅
        b1, ∅ --> X₁
        b2, ∅ --> X₂
        r, X₁ + X₂ --> 2X₂
        r, X₂ + X₁ --> 2X₁
    end
    parameters = [:r => r, :b1 => b1, :b2 => b2, :k1 => k1, :k2 => k2]
    return system, parameters
end
