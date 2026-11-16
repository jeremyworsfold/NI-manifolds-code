using NIManifoldsCode
using Random: Xoshiro, randexp

const AntsModel = NIManifoldsCode.Ants
const SOURCES = (:A, :B)
const ANTS_DENSITY_HEADROOM = 1.15

@option struct AntsModelConfig
    foragers::Int
    omega0::Float64
    omega1::Float64
    delta::Float64
end

@option struct AntsSettings
    food_steps::Int
    epsilon_panel::Float64
    epsilon_sweep::Vector{Float64}
    survival_points::Int
    survival_span::Float64
    survival_substeps::Int
    merged_horizon::Float64
    merged_seed::Int
    residence_epsilons::Vector{Float64}
    sweep_food_steps::Int
    trajectory_seed::Int
    trajectory_duration::Float64
    trajectory_step::Float64
end

@option struct AntsConfig
    model::AntsModelConfig
    quick::AntsSettings
    full::AntsSettings
end

function ants_parameters(model::AntsModelConfig, epsilon)
    return AntsModel.ModelParameters(;
        model.foragers, epsilon, model.delta, model.omega0, model.omega1
    )
end

function check_config(config::AntsConfig)
    for mode in MODES
        settings = getfield(config, mode)
        table = "[$mode] of config/ants.toml"
        sweep = settings.epsilon_sweep
        require(
            !isempty(sweep) && all(>(0), sweep) && issorted(sweep),
            "epsilon_sweep must be positive and increasing in $table",
        )
        require(settings.epsilon_panel in sweep, "epsilon_panel must be in epsilon_sweep in $table")
        require(
            settings.food_steps > 0 && settings.sweep_food_steps > 0,
            "food_steps and sweep_food_steps must be positive in $table",
        )
        require(
            settings.survival_points >= 2 &&
                settings.survival_span > 0 &&
                settings.survival_substeps > 0,
            "survival_points >= 2, survival_span > 0 and survival_substeps > 0 in $table",
        )
        require(settings.merged_horizon > 0, "merged_horizon must be positive in $table")
        require(
            settings.trajectory_duration > 0 && settings.trajectory_step > 0,
            "trajectory_duration and trajectory_step must be positive in $table",
        )
        require(
            !isempty(settings.residence_epsilons) &&
                all(e -> 0 < e <= settings.epsilon_panel, settings.residence_epsilons),
            "residence_epsilons must lie in (0, epsilon_panel] in $table",
        )
        ants_parameters(config.model, settings.epsilon_panel)
    end
end

function ant_trajectory(model, settings)
    rng = Xoshiro(settings.trajectory_seed)
    parameters = ants_parameters(model, settings.epsilon_panel)
    N, epsilon, delta, M =
        parameters.foragers, parameters.epsilon, parameters.delta, settings.food_steps
    duration, step = settings.trajectory_duration, settings.trajectory_step
    beta, end_time = 1 / epsilon, duration / delta
    t, n, food_units = 0.0, N, 0
    at_A = true
    switches = [(time=0.0, at_A=at_A)]
    grid = collect(0:step:duration)
    fraction = Vector{Float64}(undef, length(grid))
    food = similar(fraction)
    sample_index = 1
    while t < end_time
        u = food_units / M
        copying = beta * n * (N - n)
        target = M * (1 - n / N)
        rates = (
            copying + (N - n) * AntsModel.leave_B(parameters, u),
            copying + n * AntsModel.leave_A(parameters, u),
            delta * max(target - food_units, 0),
            delta * max(food_units - target, 0),
        )
        total = sum(rates)
        next_time = t + randexp(rng) / total
        while sample_index <= length(grid) && grid[sample_index] <= delta * next_time
            fraction[sample_index] = n / N
            food[sample_index] = u
            sample_index += 1
        end
        next_time > end_time && break
        t = next_time
        draw = rand(rng) * total
        if draw < rates[1]
            n += 1
        elseif draw < rates[1] + rates[2]
            n -= 1
        elseif draw < rates[1] + rates[2] + rates[3]
            food_units += 1
        else
            food_units -= 1
        end
        if draw < rates[1] + rates[2] && n in (0, N) && (n == N) != at_A
            at_A = n == N
            push!(switches, (time=delta * t, at_A=at_A))
        end
    end

    reduced_t, reduced_u = Float64[], Float64[]
    u = 0.0
    bounds = [getfield.(switches[2:end], :time); duration]
    for (visit, stop) in zip(switches, bounds)
        count = max(1, ceil(Int, (stop - visit.time) / step))
        target = visit.at_A ? 0.0 : 1.0
        for i in 0:count
            time = visit.time + (stop - visit.time) * i / count
            push!(reduced_t, time)
            push!(reduced_u, target + (u - target) * exp(-(time - visit.time)))
        end
        u = target + (u - target) * exp(-(stop - visit.time))
    end
    return (; grid, fraction, food, switches, reduced_t, reduced_u)
end

function write_ants_trajectory(output, model, settings)
    path = ant_trajectory(model, settings)
    data = joinpath(output, "data", "ants-trajectory")
    write_columns(joinpath(data, "fraction.csv"), ("t", "y"), path.grid, path.fraction)
    write_columns(joinpath(data, "food.csv"), ("t", "u"), path.grid, path.food)
    switch_times = getfield.(path.switches, :time)
    switch_states = getfield.(path.switches, :at_A)
    write_columns(
        joinpath(data, "shade.csv"),
        ("t", "s"),
        [switch_times; settings.trajectory_duration],
        [switch_states; last(switch_states)],
    )
    return write_columns(
        joinpath(data, "reduced.csv"),
        ("t", "u"),
        path.reduced_t,
        path.reduced_u,
    )
end

Base.@kwdef struct SourceResidence
    survival::Comparison{Vector{Float64}}
    mean_exact::Float64
    mean_corrected::Float64
    mean_merged::ByOrder{Float64}
    cv2_exact::Float64
    cv2_merged::ByOrder{Float64}
    survival_error::ByOrder{Float64}
end

struct FoodDensity
    lattice::Vector{Float64}
    bar_edges::Vector{Float64}
    density::Comparison{Vector{Float64}}
end

struct ResidenceCheck
    epsilon::Float64
    mean_exact::NTuple{2,Float64}
    mean_merged::NTuple{2,ByOrder{Float64}}
end

Base.@kwdef struct AntsResults
    epsilon_sweep::Vector{Float64}
    w1::ByOrder{Vector{Float64}}
    epsilon_panel::Float64
    food_steps::Int
    food::FoodDensity
    support_corrected::Vector{Float64}
    times::Vector{Float64}
    A::SourceResidence
    B::SourceResidence
    residence_checks::Vector{ResidenceCheck}
end

source_results(r::AntsResults) = (r.A, r.B)

function density_errors(model, settings)
    M = settings.sweep_food_steps
    leading = AntsModel.stationary_solution(ants_parameters(model, 0.0))
    exact = map(settings.epsilon_sweep) do epsilon
        AntsModel.exact_chain(ants_parameters(model, epsilon); food_steps=M)
    end
    corrected = map(settings.epsilon_sweep) do epsilon
        AntsModel.stationary_solution(ants_parameters(model, epsilon); correction=:local)
    end
    return ByOrder(
        [AntsModel.wasserstein_distance(leading, chain) for chain in exact],
        [AntsModel.wasserstein_distance(solution, chain) for (solution, chain) in zip(corrected, exact)],
    )
end

function food_density(settings, leading, corrected, chain)
    M = settings.food_steps
    comparison = Comparison(
        AntsModel.food_masses(chain) .* M,
        AntsModel.lattice_profile(leading, M) .* M,
        AntsModel.lattice_profile(corrected, M) .* M,
    )
    return FoodDensity(
        collect((0:M) ./ M), collect(((0:(M + 1)) .- 0.5) ./ M), comparison
    )
end

function merged_moments(parameters, settings; horizon=settings.merged_horizon)
    samples = map((:leading, :local)) do correction
        AntsModel.merged_residences(
            parameters; correction, horizon, seed=settings.merged_seed
        )
    end
    return ByOrder(samples...)
end

function source_residences(model, settings, corrected, chain)
    exact = Dict(source => AntsModel.exact_residence(chain, source) for source in SOURCES)
    times = range(
        0,
        settings.survival_span * maximum(result.mean for result in values(exact));
        length=settings.survival_points,
    )
    merged = merged_moments(
        ants_parameters(model, settings.epsilon_panel), settings; horizon=settings.merged_horizon
    )
    function source_result(source)
        sample(order) = getfield(getfield(merged, order), source)
        moments = ByOrder(
            AntsModel.sample_moments(sample(:leading)),
            AntsModel.sample_moments(sample(:corrected)),
        )
        survival = Comparison(
            AntsModel.exact_residence(
                chain, source; times, substeps=settings.survival_substeps
            ).survival,
            AntsModel.empirical_survival(sample(:leading), times),
            AntsModel.empirical_survival(sample(:corrected), times),
        )
        return SourceResidence(;
            survival,
            mean_exact=exact[source].mean,
            mean_corrected=AntsModel.residence_moments(corrected, source).mean,
            mean_merged=map(moment -> moment.mean, moments),
            cv2_exact=exact[source].cv2,
            cv2_merged=map(moment -> moment.cv2, moments),
            survival_error=errors(survival, (theory, truth) -> maximum(abs, theory .- truth)),
        )
    end
    return collect(times), source_result(:A), source_result(:B)
end

function residence_checks(model, settings, panel::AntsResults)
    return map(settings.residence_epsilons) do epsilon
        if epsilon == settings.epsilon_panel
            return ResidenceCheck(
                epsilon,
                (panel.A.mean_exact, panel.B.mean_exact),
                (panel.A.mean_merged, panel.B.mean_merged),
            )
        end
        parameters = ants_parameters(model, epsilon)
        chain = AntsModel.exact_chain(parameters; food_steps=settings.food_steps)
        merged = merged_moments(parameters, settings)
        mean_exact = Tuple(AntsModel.exact_residence(chain, source).mean for source in SOURCES)
        mean_merged = Tuple(
            ByOrder(
                AntsModel.sample_moments(getfield(merged.leading, source)).mean,
                AntsModel.sample_moments(getfield(merged.corrected, source)).mean,
            ) for source in SOURCES
        )
        return ResidenceCheck(epsilon, mean_exact, mean_merged)
    end
end

function compute_ants(model::AntsModelConfig, settings::AntsSettings)
    leading = AntsModel.stationary_solution(ants_parameters(model, 0.0))
    corrected = AntsModel.stationary_solution(
        ants_parameters(model, settings.epsilon_panel); correction=:local
    )
    chain = AntsModel.exact_chain(
        ants_parameters(model, settings.epsilon_panel); food_steps=settings.food_steps
    )
    times, source_A, source_B = source_residences(model, settings, corrected, chain)
    partial = AntsResults(;
        settings.epsilon_sweep,
        w1=density_errors(model, settings),
        settings.epsilon_panel,
        settings.food_steps,
        food=food_density(settings, leading, corrected, chain),
        support_corrected=collect(AntsModel.support(corrected)),
        times,
        A=source_A,
        B=source_B,
        residence_checks=ResidenceCheck[],
    )
    return AntsResults(;
        partial.epsilon_sweep,
        partial.w1,
        partial.epsilon_panel,
        partial.food_steps,
        partial.food,
        partial.support_corrected,
        partial.times,
        partial.A,
        partial.B,
        residence_checks=residence_checks(model, settings, partial),
    )
end

function check_ants(r::AntsResults)
    check_improves(r.w1, "the ant food-level error at the largest epsilon")
    for (source, result) in zip(SOURCES, source_results(r))
        check_improves(result.survival_error, "source $source survival")
        abs(result.mean_merged.corrected / result.mean_exact - 1) < 0.02 ||
            error("merged source $source residence mean differs from exact by more than 2%")
    end
    return nothing
end

function write_ants_data(output, settings, r::AntsResults)
    (; lattice, bar_edges, density) = r.food
    data = joinpath(output, "data", "ants-correction")
    write_columns(
        joinpath(data, "errors-epsilon.csv"),
        ("epsilon", "leading", "corrected"),
        r.epsilon_sweep,
        r.w1.leading,
        r.w1.corrected,
    )
    write_columns(
        joinpath(data, "density.csv"),
        ("u", "leading", "corrected"),
        lattice,
        density.leading,
        density.corrected,
    )
    write_columns(
        joinpath(data, "exact.csv"), ("u", "density"), bar_edges[1:(end - 1)], density.exact
    )
    write_columns(
        joinpath(data, "survival.csv"),
        ("t", "exactA", "leadingA", "correctedA", "exactB", "leadingB", "correctedB"),
        r.times,
        r.A.survival.exact,
        r.A.survival.leading,
        r.A.survival.corrected,
        r.B.survival.exact,
        r.B.survival.leading,
        r.B.survival.corrected,
    )
    return write_tex_macros(
        joinpath(data, "params.tex"),
        (
            "epsPanel" => r.epsilon_panel,
            "foodSteps" => r.food_steps,
            "sweepFoodSteps" => settings.sweep_food_steps,
            "densityMax" => ANTS_DENSITY_HEADROOM * maximum(density.exact),
            "survivalSpan" => last(r.times),
            "meanResidenceA" => r.A.mean_exact,
            "meanResidenceB" => r.B.mean_exact,
        ),
    )
end

function ants_metrics(r::AntsResults)
    function source_table(result)
        return Dict{String,Any}(
            "mean_exact" => result.mean_exact,
            "mean_corrected" => result.mean_corrected,
            "mean_merged" => metrics_table(result.mean_merged),
            "cv2_exact" => result.cv2_exact,
            "cv2_merged" => metrics_table(result.cv2_merged),
            "survival_error" => metrics_table(result.survival_error),
        )
    end
    function check_table(check)
        table = Dict{String,Any}("epsilon" => check.epsilon)
        for (i, source) in enumerate(SOURCES)
            exact, merged = check.mean_exact[i], check.mean_merged[i]
            table[String(source)] = Dict{String,Any}(
                "mean_exact" => exact,
                "mean_merged" => metrics_table(merged),
                "relative_error" => metrics_table(map(value -> value / exact - 1, merged)),
            )
        end
        return table
    end
    return Dict{String,Any}(
        "w1" => metrics_table(r.w1),
        "panel" => Dict{String,Any}(
            "support_corrected" => r.support_corrected,
            "A" => source_table(r.A),
            "B" => source_table(r.B),
        ),
        "residence_check" => Dict(
            epsilon_key(check.epsilon) => check_table(check) for check in r.residence_checks
        ),
    )
end

function reproduce_ants(mode::Symbol)
    config = load_config(AntsConfig, "ants")
    model, settings = config.model, getfield(config, mode)
    output = output_dir("ants", mode)
    results = compute_ants(model, settings)
    check_ants(results)
    write_ants_data(output, settings, results)
    write_ants_trajectory(output, model, settings)
    save_metrics(output, mode, model, settings, ants_metrics(results))
    return output
end
