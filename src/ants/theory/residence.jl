# Residence times in the reduced PDMP. A visit starts from the stationary flux into a source
# and ends at the first event of the switching clock along the deterministic food-level flow.

function _drift(p::ModelParameters, c, u, source::Symbol)
    correction = c * p.epsilon * p.delta * switching_rate(p, u, source) / 2
    return source === :A ? -p.delta * u + correction : p.delta * (1 - u) - correction
end

"""Entry weights and food levels for visits to `source`."""
function _visit(s::StationarySolution, source::Symbol)
    _check_source(source)
    p = s.parameters
    shift = _order(s.correction) * p.epsilon * p.delta / 2
    x, m = s.nodes, s.masses
    weights, entries = if source === :A
        m .* x .* leave_B.(Ref(p), x), x .+ shift
    else
        m .* (1 .- x) .* leave_A.(Ref(p), x), x .- shift
    end
    return weights ./ sum(weights), entries
end

function _visit_problem(s::StationarySolution, source::Symbol, horizon)
    weights, entries = _visit(s, source)
    p, c, n = s.parameters, _order(s.correction), length(entries)
    function rates!(dy, y, _, t)
        for i in 1:n
            u, G = y[i], exp(-y[n + i])
            dy[i] = _drift(p, c, u, source)
            dy[n + i] = switching_rate(p, u, source)
            dy[2n + i] = G
            dy[3n + i] = t * G
        end
    end
    return weights, ODEProblem(rates!, [entries; zeros(3n)], (0.0, float(horizon)))
end

"""Survival function of a complete visit to `source` at the sorted `times`."""
function survival(s::StationarySolution, source::Symbol, times; reltol=1e-10, abstol=1e-12)
    issorted(times) && first(times) >= 0 || throw(ArgumentError("times must be sorted and >= 0"))
    weights, problem = _visit_problem(s, source, last(times))
    n = length(weights)
    path = solve(problem, Tsit5(); saveat=times, save_idxs=(n + 1):(2n), reltol, abstol)
    return [dot(weights, exp.(-integrated)) for integrated in path.u]
end

"""Mean, second moment and squared coefficient of variation of visits to `source`."""
function residence_moments(
    s::StationarySolution, source::Symbol; reltol=1e-10, abstol=1e-12, tolerance=1e-14
)
    weights, problem = _visit_problem(s, source, Inf)
    n = length(weights)
    finished(y, _, _) = all(>(-log(tolerance)), view(y, (n + 1):(2n)))
    path = solve(
        problem,
        Tsit5();
        callback=DiscreteCallback(finished, terminate!),
        save_everystep=false,
        reltol,
        abstol,
    )
    y = last(path.u)
    mean_value = dot(weights, view(y, (2n + 1):(3n)))
    second = 2 * dot(weights, view(y, (3n + 1):(4n)))
    return (; mean=mean_value, second, cv2=second / mean_value^2 - 1)
end
