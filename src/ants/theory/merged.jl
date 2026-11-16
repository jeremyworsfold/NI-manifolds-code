# The exact chain counts a visit from first arrival at one unanimous source to first arrival
# at the other. A reduced visit shorter than the fast recruitment transit is therefore merged
# with its neighbours before residence statistics are compared.

"""Duration of a successful neutral-recruitment excursion conditioned on fixation."""
function transit_time(rng::AbstractRNG, N::Integer, beta::Real)
    n, t = 1, 0.0
    while n < N
        t += randexp(rng) / (2beta * n * (N - n))
        n += rand(rng) < (n + 1) / (2n) ? 1 : -1
    end
    return t
end

"""Switch times and the source entered (`true` means A) along a reduced path."""
function _switches(p::ModelParameters, c, horizon, rng)
    shift = c * p.epsilon * p.delta / 2
    at_A = Ref(true)
    source() = at_A[] ? :A : :B
    times, entered_A = Float64[], Bool[]
    flow!(du, u, _, _) = (du[1] = _drift(p, c, u[1], source()))
    rate(u, _, _) = switching_rate(p, u[1], source())
    function switch!(integrator)
        integrator.u[1] += at_A[] ? -shift : shift
        at_A[] = !at_A[]
        push!(times, integrator.t)
        return push!(entered_A, at_A[])
    end
    problem = JumpProblem(
        ODEProblem(flow!, [0.0], (0.0, float(horizon))),
        VariableRateJump(rate, switch!);
        vr_aggregator=VR_Direct(),
        rng,
    )
    path = solve(
        problem, Tsit5(); save_everystep=false, reltol=1e-8, abstol=1e-10, maxiters=typemax(Int)
    )
    last(path.t) == horizon || error("reduced PDMP path stopped at t = $(last(path.t))")
    return times, entered_A
end

function merged_residences(
    p::ModelParameters;
    correction::Symbol=:leading,
    horizon::Real=1e6,
    seed::Integer=1,
    burn::Real=50.0,
)
    path_rng, transit_rng = Xoshiro(seed), Xoshiro(seed + 1)
    times, entered_A = _switches(p, _order(correction), horizon, path_rng)
    length(times) >= 2 || error("reduced path produced too few switches")
    beta = 1 / p.epsilon
    kept_times, kept_sources = [times[1]], [entered_A[1]]
    i = 2
    while i < length(times)
        if p.epsilon > 0 &&
           times[i + 1] - times[i] < transit_time(transit_rng, p.foragers, beta)
            i += 2
        else
            push!(kept_times, times[i])
            push!(kept_sources, entered_A[i])
            i += 1
        end
    end
    source_A, source_B = Float64[], Float64[]
    for k in 1:(length(kept_times) - 1)
        kept_times[k] > burn || continue
        push!(kept_sources[k] ? source_A : source_B, kept_times[k + 1] - kept_times[k])
    end
    return (; A=source_A, B=source_B)
end

function sample_moments(residences::AbstractVector)
    mean_value, second = mean(residences), mean(abs2, residences)
    return (; mean=mean_value, second, cv2=second / mean_value^2 - 1)
end

empirical_survival(residences::AbstractVector, times) = 1 .- ecdf(residences).(times)
