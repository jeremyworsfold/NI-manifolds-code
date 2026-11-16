@testset "Ant-foraging parameters and rates" begin
    p = AntsModel.ModelParameters(
        foragers=20, epsilon=0.1, delta=0.22, omega0=0.02, omega1=0.2
    )
    @test AntsModel.leave_A(p, 1.0) == p.omega0
    @test AntsModel.leave_A(p, 0.0) == p.omega0 + p.omega1
    @test AntsModel.leave_B(p, 0.0) == p.omega0
    @test AntsModel.leave_B(p, 1.0) == p.omega0 + p.omega1
    @test AntsModel.switching_rate(p, 0.3, :A) == AntsModel.leave_A(p, 0.3)
    @test AntsModel.switching_rate(p, 0.3, :B) == AntsModel.leave_B(p, 0.3)
    @test_throws ArgumentError AntsModel.switching_rate(p, 0.3, :plus)
    @test_throws ArgumentError AntsModel.ModelParameters(foragers=1)
    @test_throws ArgumentError AntsModel.ModelParameters(epsilon=-0.1)
    @test_throws ArgumentError AntsModel.ModelParameters(delta=0.0)
    @test_throws ArgumentError AntsModel.ModelParameters(omega0=0.0)
    @test_throws ArgumentError AntsModel.ModelParameters(omega1=-0.1)
end

@testset "Ant-foraging stationary theory" begin
    p = AntsModel.ModelParameters(
        foragers=20, epsilon=0.1, delta=0.22, omega0=0.02, omega1=0.2
    )
    leading = AntsModel.stationary_solution(
        AntsModel.ModelParameters(
            foragers=p.foragers,
            epsilon=0.0,
            delta=p.delta,
            omega0=p.omega0,
            omega1=p.omega1,
        )
    )
    corrected = AntsModel.stationary_solution(p; correction=:local)
    @test AntsModel.support(leading) == (0.0, 1.0)
    lower, upper = AntsModel.support(corrected)
    @test lower ≈ 1 - upper rtol = 1e-12
    @test 0 < lower < 0.5 < upper < 1
    @test collect(AntsModel.source_probabilities(leading)) ≈ [0.5, 0.5] rtol = 1e-12
    @test collect(AntsModel.source_probabilities(corrected)) ≈ [0.5, 0.5] rtol = 1e-12
    @test AntsModel.density(leading, 0.25) ≈ AntsModel.density(leading, 0.75) rtol = 1e-12
    @test AntsModel.source_density(leading, 0.3, :A) ≈
        AntsModel.source_density(leading, 0.7, :B) rtol = 1e-12
    @test AntsModel.cumulative(leading, 0.5) ≈ 0.5 rtol = 1e-12
    @test sum(AntsModel.lattice_profile(corrected, 80)) ≈ 1 rtol = 1e-11
    @test_throws ArgumentError AntsModel.stationary_solution(p; correction=:reset)
    @test_throws ArgumentError AntsModel.source_density(leading, 0.5, :plus)
end

@testset "Ant-foraging exact chain" begin
    p = AntsModel.ModelParameters(
        foragers=10, epsilon=0.1, delta=0.22, omega0=0.02, omega1=0.2
    )
    chain = AntsModel.exact_chain(p; food_steps=200)
    @test norm(chain.generator * ones(length(chain.stationary)), Inf) < 1e-10
    @test norm(transpose(chain.generator) * chain.stationary, Inf) < 1e-10
    @test sum(AntsModel.food_masses(chain)) ≈ 1
    @test collect(AntsModel.source_probabilities(chain)) ≈ [0.5, 0.5] rtol = 1e-10
    exact_A = AntsModel.exact_residence(chain, :A)
    exact_B = AntsModel.exact_residence(chain, :B)
    @test exact_A.mean ≈ exact_B.mean rtol = 1e-9
    @test exact_A.cv2 ≈ exact_B.cv2 rtol = 1e-8

    leading = AntsModel.stationary_solution(
        AntsModel.ModelParameters(
            foragers=p.foragers,
            epsilon=0.0,
            delta=p.delta,
            omega0=p.omega0,
            omega1=p.omega1,
        )
    )
    corrected = AntsModel.stationary_solution(p; correction=:local)
    @test AntsModel.wasserstein_distance(corrected, chain) <
        AntsModel.wasserstein_distance(leading, chain)
    @test_throws ArgumentError AntsModel.exact_chain(p; food_steps=0)
end

@testset "Ant-foraging residence theory" begin
    p = AntsModel.ModelParameters(
        foragers=10, epsilon=0.05, delta=0.22, omega0=0.02, omega1=0.2
    )
    corrected = AntsModel.stationary_solution(p; correction=:local, nodes=100)
    moments_A = AntsModel.residence_moments(corrected, :A)
    moments_B = AntsModel.residence_moments(corrected, :B)
    @test moments_A.mean ≈ moments_B.mean rtol = 1e-10
    @test moments_A.cv2 ≈ moments_B.cv2 rtol = 1e-10
    times = range(0, 20; length=41)
    @test AntsModel.survival(corrected, :A, times) ≈
        AntsModel.survival(corrected, :B, times) rtol = 1e-10
end
