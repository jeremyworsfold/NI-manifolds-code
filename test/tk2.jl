@testset "TK2 microscopic process" begin
    reduced = TK2Model.ReducedParameters(kappa=0.6, beta=0.0, epsilon=0.1)
    model = TK2Model.realize(reduced; population=20)
    @test TK2Model.reduce(model) == reduced

    chain = TK2Model.exact_chain(model; cutoff=65)
    @test norm(chain.generator * ones(length(chain.stationary)), Inf) < 1e-10
    @test norm(transpose(chain.generator) * chain.stationary, Inf) < 1e-11
    @test TK2Model.boundary_mass(chain) < 1e-3
    @test sum(TK2Model.population_masses(chain)) ≈ 1
    @test sum(TK2Model.branch_probabilities(chain)) ≈ 1

    # The default cutoff is `sigma` Poisson widths beyond N/(1 - |κ|), plus `extra`.
    @test TK2Model.default_cutoff(model; sigma=8, extra=20) ==
        ceil(Int, 20 / 0.4 + 8 * sqrt(20 / 0.4)) + 20
end

function beta_pushforward_moment(kappa, beta_bias, order)
    cplus, cminus = (1 + beta_bias) / 2, (1 - beta_bias) / 2
    lower, upper = 1 / (1 + kappa), 1 / (1 - kappa)
    integrate(f) = solve(IntegralProblem((x, _) -> f(x), (0.0, 1.0)), QuadGKJL(); reltol=1e-11).u
    left = integrate() do v
        u = 0.5 * v^(1 / cminus)
        (1 - u)^(cplus - 1) * (lower + (upper - lower) * u)^order
    end
    right = integrate() do v
        u = 1 - 0.5 * v^(1 / cplus)
        u^(cminus - 1) * (lower + (upper - lower) * u)^order
    end
    return (0.5^cminus / cminus * left + 0.5^cplus / cplus * right) / exp(logbeta(cminus, cplus))
end

@testset "TK2 stationary theory" begin
    for (kappa, beta, epsilon) in ((0.4, -0.2, 0.2), (0.6, 0.0, 0.1))
        parameters = TK2Model.ReducedParameters(; kappa, beta, epsilon)
        leading = TK2Model.stationary_solution(parameters; correction=:leading)
        moments = TK2Model.stationary_moments(leading)
        @test moments.mean ≈ beta_pushforward_moment(kappa, beta, 1) rtol = 1e-9
        @test moments.second ≈ beta_pushforward_moment(kappa, beta, 2) rtol = 1e-9

        local_solution = TK2Model.stationary_solution(parameters; correction=:local)
        cplus, cminus = (1 + beta) / 2, (1 - beta) / 2
        plus_flux =
            cminus * (1 + epsilon * kappa) * TK2Model._branch_integral(local_solution, :plus, -1)
        minus_flux =
            cplus * (1 - epsilon * kappa) * TK2Model._branch_integral(local_solution, :minus, -1)
        @test plus_flux ≈ minus_flux rtol = 1e-12
    end

    zero = TK2Model.ReducedParameters(kappa=0.4, beta=-0.2, epsilon=0.0)
    leading = TK2Model.stationary_solution(zero; correction=:leading)
    for correction in (:local, :reset)
        solution = TK2Model.stationary_solution(zero; correction)
        @test collect(TK2Model.branch_probabilities(solution)) ≈
            collect(TK2Model.branch_probabilities(leading)) rtol = 1e-12
        @test TK2Model.switching_flux(solution) ≈ TK2Model.switching_flux(leading) rtol = 1e-12
    end
    flat = TK2Model.stationary_solution(
        TK2Model.ReducedParameters(kappa=0.0, beta=0.2, epsilon=0.2); correction=:reset
    )
    @test collect(TK2Model.branch_probabilities(flat)) ≈ [0.6, 0.4]
    @test TK2Model.stationary_moments(flat) == (mean=1.0, second=1.0, variance=0.0)

    # Near κ = 0 the corrected formulas would divide by κ, so the leading law is used.
    tiny = TK2Model.ReducedParameters(kappa=1e-9, beta=0.2, epsilon=0.1)
    @test collect(
        TK2Model.branch_probabilities(TK2Model.stationary_solution(tiny; correction=:reset))
    ) ≈ [0.6, 0.4] rtol = 1e-6
    @test_throws ArgumentError TK2Model.ReducedParameters(kappa=-0.4)
end

@testset "TK2 microscopic correction validation" begin
    microscopic = TK2Model.realize(
        TK2Model.ReducedParameters(kappa=0.6, beta=0.0, epsilon=0.05); population=15
    )
    exact_plus = first(TK2Model.branch_probabilities(TK2Model.exact_chain(microscopic; cutoff=75)))
    leading_plus = first(
        TK2Model.branch_probabilities(
            TK2Model.stationary_solution(TK2Model.reduce(microscopic); correction=:leading)
        ),
    )
    corrected_plus = first(
        TK2Model.branch_probabilities(
            TK2Model.stationary_solution(TK2Model.reduce(microscopic); correction=:reset)
        ),
    )
    @test abs(corrected_plus - exact_plus) < abs(leading_plus - exact_plus) / 10
end

@testset "TK2 stationary densities" begin
    parameters = TK2Model.ReducedParameters(kappa=0.4, beta=-0.2, epsilon=0.2)
    for (correction, type) in (
        (:leading, TK2Model.LeadingPDMP), (:local, TK2Model.LocalPDMP), (:reset, TK2Model.ResetPDMP)
    )
        solution = TK2Model.stationary_solution(parameters; correction)
        @test solution isa type
        @test solution.parameters == parameters
        lower, upper = TK2Model.support(solution)
        # The densities have integrable power-law singularities at both ends of the support,
        # so stop `INSET` short of them, which loses a relative mass of order 1e-6.
        integrate(f) = solve(
            IntegralProblem((x, _) -> f(x), (lower + TK2Model.INSET, upper - TK2Model.INSET)),
            QuadGKJL();
            reltol=1e-10,
        ).u
        plus, minus = TK2Model.branch_probabilities(solution)
        @test integrate(x -> TK2Model.branch_density(solution, x, :plus)) ≈ plus rtol = 1e-4
        @test integrate(x -> TK2Model.branch_density(solution, x, :minus)) ≈ minus rtol = 1e-4
        @test integrate(x -> TK2Model.density(solution, x)) ≈ 1 rtol = 1e-4
        @test_throws ArgumentError TK2Model.branch_density(solution, 1.0, :left)
    end
    # At κ = 0 every correction is the leading law.
    flat = TK2Model.ReducedParameters(kappa=0.0, beta=0.0, epsilon=0.1)
    @test TK2Model.stationary_solution(flat; correction=:reset) isa TK2Model.LeadingPDMP
end

@testset "TK2 figure helpers" begin
    parameters = TK2Model.ReducedParameters(kappa=0.4, beta=-0.2, epsilon=0.0)
    leading = TK2Model.stationary_solution(parameters; correction=:leading)
    plus, minus = TK2Model.branch_probabilities(leading)
    @test TK2Model.mean_branch_coordinate(leading) == plus - minus
    flat = TK2Model.stationary_solution(
        TK2Model.ReducedParameters(kappa=0.0, beta=0.2, epsilon=0.2); correction=:reset
    )
    @test TK2Model.mean_branch_coordinate(flat) ≈ 0.2

    # The L1 distance is zero against the smoothed density's own lattice masses, and grows by
    # exactly the size of any perturbation of them.
    V = 40
    grid = (0:80) ./ V
    masses = TK2Model.finite_size_profile(leading, V; grid).total ./ V
    @test TK2Model.density_l1(leading, masses, V) ≈ 0 atol = 1e-15
    shift = [isodd(n) ? 1e-4 : -1e-4 for n in eachindex(masses)]
    @test TK2Model.density_l1(leading, masses .+ shift, V) ≈ sum(abs, shift) rtol = 1e-9

    path = TK2Model.sample_path(leading, 1.0, 5; rng=Xoshiro(3))
    @test length(path) == 5
    @test [segment.branch for segment in path] == [:plus, :minus, :plus, :minus, :plus]
    lower, upper = TK2Model.support(leading)
    for (segment, next) in zip(path[1:(end - 1)], path[2:end])
        @test last(segment.t) ≈ first(next.t)  # residences tile the time axis
        @test last(segment.x) ≈ first(next.x)  # the state is continuous at a switch
    end
    for segment in path
        target = segment.branch === :plus ? lower : upper  # the fixed point 1/k of the branch
        @test issorted(abs.(segment.x .- target); rev=true)
    end
    @test [last(s.x) for s in path] == [last(s.x) for s in TK2Model.sample_path(leading, 1.0, 5; rng=Xoshiro(3))]
end
