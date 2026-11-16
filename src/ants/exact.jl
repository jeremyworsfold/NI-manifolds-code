# Exact labelled chain on (n, m, ell): n foragers are at A, m of M food units remain there,
# and ell records the unanimous source visited last. Neutral recruitment has beta = 1/epsilon.
# Food moves by a bounded upwind discretisation toward M(1-n/N), so it is consumed at A and
# recovers at B without leaving [0, M].

struct ExactChain
    parameters::ModelParameters
    food_steps::Int
    generator::SparseMatrixCSC{Float64,Int}
    stationary::Vector{Float64}
end

function _unlabelled_generator(p::ModelParameters, food_steps::Integer)
    N, beta, M, delta = p.foragers, 1 / p.epsilon, food_steps, p.delta
    states = LinearIndices((N + 1, M + 1))
    rows, columns, rates = Int[], Int[], Float64[]
    function move!(i, j, rate)
        rate > 0 || return nothing
        push!(rows, i, i)
        push!(columns, j, i)
        push!(rates, rate, -rate)
        return nothing
    end
    for m in 0:M, n in 0:N
        i, u = states[n + 1, m + 1], m / M
        copying = beta * n * (N - n)
        n < N && move!(i, states[n + 2, m + 1], copying + (N - n) * leave_B(p, u))
        n > 0 && move!(i, states[n, m + 1], copying + n * leave_A(p, u))
        target = M * (1 - n / N)
        m < M && move!(i, states[n + 1, m + 2], delta * max(target - m, 0))
        m > 0 && move!(i, states[n + 1, m], delta * max(m - target, 0))
    end
    return sparse(rows, columns, rates, length(states), length(states))
end

_state_index(N, food_steps) = LinearIndices((food_steps + 1, N + 1, 2))

function exact_chain(p::ModelParameters; food_steps::Integer=500)
    p.epsilon > 0 || throw(ArgumentError("the exact chain needs epsilon > 0"))
    food_steps > 0 || throw(ArgumentError("food_steps must be positive"))
    N, M = p.foragers, Int(food_steps)
    base = _unlabelled_generator(p, M)
    base_states = CartesianIndices((N + 1, M + 1))
    index = _state_index(N, M)
    rows, columns, rates = Int[], Int[], Float64[]
    for (i, j, rate) in zip(findnz(base)...)
        n, food = Tuple(base_states[i]) .- 1
        next_n, next_food = Tuple(base_states[j]) .- 1
        for label in 0:1
            next_label = if next_n == n
                label
            elseif next_n == 0
                0
            elseif next_n == N
                1
            else
                label
            end
            push!(rows, index[food + 1, n + 1, label + 1])
            push!(columns, index[next_food + 1, next_n + 1, next_label + 1])
            push!(rates, rate)
        end
    end
    generator = sparse(rows, columns, rates, length(index), length(index))
    return ExactChain(p, M, generator, stationary_vector(generator))
end

"""Stationary masses at food levels m/M, m = 0, ..., M."""
food_masses(chain::ExactChain) =
    vec(sum(reshape(chain.stationary, chain.food_steps + 1, :); dims=2))

function wasserstein_distance(s::StationarySolution, chain::ExactChain)
    return wasserstein_distance(s, food_masses(chain), chain.food_steps)
end

function _source_states(chain::ExactChain, source::Symbol)
    _check_source(source)
    index = _state_index(chain.parameters.foragers, chain.food_steps)
    return vec(index[:, :, source === :A ? 2 : 1])
end

"""Stationary exact-chain occupancies `(P_A, P_B)`."""
function source_probabilities(chain::ExactChain)
    return Tuple(sum(chain.stationary[_source_states(chain, source)]) for source in (:A, :B))
end

function exact_residence(chain::ExactChain, source::Symbol; times=nothing, substeps=4)
    _check_source(source)
    N, M, Q = chain.parameters.foragers, chain.food_steps, chain.generator
    start, previous, label = source === :A ? (N, N - 1, 0) : (0, 1, 1)
    transient = _source_states(chain, source)
    entry = zeros(size(Q, 1))
    index = _state_index(N, M)
    for food in 0:M
        i = index[food + 1, previous + 1, label + 1]
        j = index[food + 1, start + 1, 2 - label]
        entry[j] = chain.stationary[i] * Q[i, j]
    end
    flux = sum(entry)
    weights = entry[transient] ./ flux
    QT = Q[transient, transient]
    F = lu(QT)
    tau1 = F \ (-ones(length(transient)))
    tau2 = F \ (-2 .* tau1)
    mean_value, second = dot(weights, tau1), dot(weights, tau2)
    result = (;
        mean=mean_value,
        second,
        cv2=second / mean_value^2 - 1,
        flux,
        occupancy=sum(chain.stationary[transient]),
    )
    isnothing(times) && return result
    return (; result..., survival=_transient_survival(QT, weights, times, substeps))
end

function _transient_survival(Q, q0, times, substeps)
    first(times) == 0 || throw(ArgumentError("times must start at zero"))
    problem = ODEProblem(MatrixOperator(sparse(transpose(Q))), q0, (0.0, last(times)))
    solver = TRBDF2(; linsolve=KLUFactorization())
    path = solve(problem, solver; adaptive=false, dt=step(times) / substeps, saveat=times)
    return sum.(path.u)
end
