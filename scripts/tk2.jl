using NIManifoldsCode
using Random

const TK2Model = NIManifoldsCode.TK2

# --- Configuration: config/tk2.toml -------------------------------------------------

"""
The `[model]` table of `config/tk2.toml`: the reduced TK2 parameters, the same in both modes.

The manuscript's TK2 figures compare the stationary density and the mean branch coordinate
y = P₊ - P₋ of the exact chain with the leading and O(ε) reset PDMP, against the timescale
separation ε.

# Fields
- `kappa`: decay bias κ, with 0 ≤ κ < 1.
- `beta`: influx bias β, with |β| ≤ 1.
"""
@option struct TK2ModelConfig
    kappa::Float64
    beta::Float64
end

"""
A `[quick]` or `[full]` table of `config/tk2.toml`: resolutions and sweeps for one mode. The
panels are those of the correction figure.

# Fields
Stationary densities, panel (b):
- `population`: volume V of the density panels and of the exact samples in (a).
- `points`: grid points of the smoothed PDMP densities.
- `epsilon_rows`: the two ε of the density rows, which must lie on the exact lattice of (a).

Mean branch coordinate, panel (a):
- `epsilon_plot_min`, `epsilon_plot_max`: the ε range.
- `curve_points`: points of the O(ε) theory curve, log-spaced.
- `exact_points_per_decade`: exact samples per decade, on a log lattice anchored at
  `epsilon_rows[1]`.

Exact-chain cutoff, passed to `TK2.default_cutoff`:
- `cutoff_sigma`: Poisson widths beyond the upper branch fixed point.
- `cutoff_extra`: extra states beyond that.

Density-error convergence, panels (c, d). The lattice L1 error against the exact chain is swept
in ε at each of two V, and in μ = 1/V at each of two small ε, where the ε error does not mask
the finite-volume one:
- `error_epsilons`: the ε of panel (c).
- `error_populations_epsilon`: its two V.
- `error_populations_mu`: the V of panel (d).
- `error_mu_epsilons`: its two ε.
"""
@option struct TK2Settings
    population::Int
    points::Int
    epsilon_plot_min::Float64
    epsilon_plot_max::Float64
    curve_points::Int
    exact_points_per_decade::Int
    epsilon_rows::Vector{Float64}
    cutoff_sigma::Float64
    cutoff_extra::Int
    error_epsilons::Vector{Float64}
    error_populations_epsilon::Vector{Int}
    error_populations_mu::Vector{Int}
    error_mu_epsilons::Vector{Float64}
end

"""`config/tk2.toml`: a shared `[model]` table and one table per mode."""
@option struct TK2Config
    model::TK2ModelConfig
    quick::TK2Settings
    full::TK2Settings
end

function check_config(config::TK2Config)
    (; kappa, beta) = config.model
    TK2Model.ReducedParameters(; kappa, beta, epsilon=0.0)  # rejects κ ∉ [0, 1) and |β| > 1
    for mode in MODES
        settings = getfield(config, mode)
        table = "[$mode] of config/tk2.toml"
        require(
            settings.population > 0 && settings.points > 1 && settings.curve_points > 1,
            "population, points and curve_points must be positive in $table",
        )
        require(
            0 < settings.epsilon_plot_min < settings.epsilon_plot_max,
            "need 0 < epsilon_plot_min < epsilon_plot_max in $table",
        )
        require(
            settings.exact_points_per_decade > 0,
            "exact_points_per_decade must be positive in $table",
        )
        # The writers index rows a/b of panel (b) and the two curves of panels (c, d).
        for field in (:epsilon_rows, :error_populations_epsilon, :error_mu_epsilons)
            require(
                length(getfield(settings, field)) == 2, "$field must have two entries in $table"
            )
        end
        tk2_exact_epsilons(settings)  # throws unless the row ε lie on the exact lattice
    end
end

# --- Results -----------------------------------------------------------------------

"""
Correction panel (a): the mean branch coordinate y = P₊ - P₋ against ε.

# Fields
- `epsilon_theory`, `corrected_y`: the O(ε) PDMP on a fine log grid.
- `epsilon_exact`, `exact_y`: the exact chain on its log lattice (`tk2_exact_epsilons`).
"""
Base.@kwdef struct MeanCurves
    epsilon_theory::Vector{Float64}
    corrected_y::Vector{Float64}
    epsilon_exact::Vector{Float64}
    exact_y::Vector{Float64}
end

"""
One row of correction panel (b): the stationary density of s = n/V at one ε.

# Fields
- `epsilon`: the row's ε.
- `grid`: the s grid of the LNA-smoothed PDMP densities.
- `density`: on `grid`, the exact lattice density sampled at the nearest lattice point, and both
  smoothed PDMP densities.
- `bar_edges`, `bar_heights`: the exact lattice density as bars, for a const plot.
- `mean_y_exact`, `mean_y_corrected`: the mean branch coordinate of the chain and the O(ε) PDMP.
- `l1`: lattice L1 errors of the densities, as in panels (c, d).
- `sup`: sup errors on `grid`, relative to the exact peak.
"""
Base.@kwdef struct DensityRow
    epsilon::Float64
    grid::Vector{Float64}
    density::Comparison{Vector{Float64}}
    bar_edges::Vector{Float64}
    bar_heights::Vector{Float64}
    mean_y_exact::Float64
    mean_y_corrected::Float64
    l1::ByOrder{Float64}
    sup::ByOrder{Float64}
end

"""Correction panel (c): Wasserstein-1 (plotted) and lattice L1 density errors against ε at one volume `population`."""
struct EpsilonSweep
    population::Int
    epsilons::Vector{Float64}
    w1::ByOrder{Vector{Float64}}
    l1::ByOrder{Vector{Float64}}
end

"""Correction panel (d): Wasserstein-1 (plotted) and lattice L1 density errors against μ = 1/V at one `epsilon`."""
struct MuSweep
    epsilon::Float64
    mus::Vector{Float64}
    w1::ByOrder{Vector{Float64}}
    l1::ByOrder{Vector{Float64}}
end

"""
Everything the two TK2 figures show.

# Fields
- `leading`: the leading-order PDMP, and `leading_y` its mean branch coordinate.
- `cutoff`: the exact-chain cutoff at the panel volume.
- `maximum_boundary_mass`: the largest stationary mass on the cutoff over every exact solve.
- `path`: the schematic path of the leading-order figure.
- `means`, `rows`, `epsilon_sweeps`, `mu_sweeps`: correction panels (a), (b), (c) and (d).
"""
Base.@kwdef struct TK2Results{L}
    leading::L
    leading_y::Float64
    cutoff::Int
    maximum_boundary_mass::Float64
    path::Vector{TK2Model.PathSegment}
    means::MeanCurves
    rows::Vector{DensityRow}
    epsilon_sweeps::Vector{EpsilonSweep}
    mu_sweeps::Vector{MuSweep}
end

# --- Theory and exact chain ---------------------------------------------------------

function tk2_leading(model::TK2ModelConfig)
    return TK2Model.stationary_solution(
        TK2Model.ReducedParameters(; model.kappa, model.beta, epsilon=0.0); correction=:leading
    )
end

function tk2_corrected(model::TK2ModelConfig, epsilon)
    return TK2Model.stationary_solution(
        TK2Model.ReducedParameters(; model.kappa, model.beta, epsilon); correction=:reset
    )
end

"""
Exact-chain cutoff at volume `population`: `cutoff_sigma` Poisson widths and `cutoff_extra`
states beyond the upper branch fixed point.
"""
function tk2_cutoff(model, settings, population=settings.population)
    return TK2Model.default_cutoff(
        TK2Model.ModelParameters(; population, decay_bias=model.kappa);
        sigma=settings.cutoff_sigma,
        extra=settings.cutoff_extra,
    )
end

"""
What the figures need from one exact-chain solve at `reduced` and volume `population`: the mean
branch coordinate, the stationary masses of the total population n = 0, 1, …, and the mass on
the cutoff, which must be negligible.
"""
function tk2_exact_moments(reduced, population::Int, cutoff::Int)
    chain = TK2Model.exact_chain(TK2Model.realize(reduced; population); cutoff)
    return (;
        mean_y=TK2Model.mean_branch_coordinate(chain),
        masses=TK2Model.population_masses(chain),
        boundary_mass=TK2Model.boundary_mass(chain),
    )
end

const TK2_EXACT_CACHE = joinpath(REPOSITORY_ROOT, "output", "cache", "tk2-exact")
# Bump when the exact chain changes, so stale cache entries are never read.
const TK2_EXACT_VERSION = 1

"""
`tk2_exact_moments`, cached on disk in `output/cache/tk2-exact/` under its parameters and
cutoff, so that redrawing the figures skips the exact solves that dominate the run time. Bump
`TK2_EXACT_VERSION` after changing the exact chain itself.
"""
function tk2_exact_moments_cached(reduced, population::Int, cutoff::Int)
    (; kappa, beta, epsilon) = reduced
    name = "fsp-rectangular_v$(TK2_EXACT_VERSION)_kappa=$(kappa)_beta=$(beta)_epsilon=$(epsilon)_V=$(population)_cutoff=$(cutoff).toml"
    path = joinpath(TK2_EXACT_CACHE, name)
    if isfile(path)
        cached = TOML.parsefile(path)
        return (;
            mean_y=cached["mean_y"],
            masses=Float64.(cached["masses"]),
            boundary_mass=cached["boundary_mass"],
        )
    end
    exact = tk2_exact_moments(reduced, population, cutoff)
    mkpath(TK2_EXACT_CACHE)
    # Write then rename, so an interrupted run never leaves a partial entry.
    save_toml("$path.partial", Dict(String(k) => v for (k, v) in pairs(exact)))
    mv("$path.partial", path; force=true)
    return exact
end

"""
Exact-chain solves on demand: `exact(ε)` at the panel volume, `exact(ε, V)` at volume `V`.
Each is memoised in `solved` (and cached on disk), so the cutoff-mass check sees all of them.
"""
function tk2_exact_solver(model, settings)
    solved = Dict{Tuple{Float64,Int},NamedTuple}()
    function exact(epsilon, V::Int=settings.population)
        get!(solved, (Float64(epsilon), V)) do
            reduced = TK2Model.ReducedParameters(; model.kappa, model.beta, epsilon)
            return tk2_exact_moments_cached(reduced, V, tk2_cutoff(model, settings, V))
        end
    end
    return exact, solved
end

# --- Leading-order figure: a schematic PDMP path ------------------------------------

"""
A short, legible path of the leading-order PDMP from s = 1: the first seed whose `nres`
residences each last between `tau_min` and `tau_max`. Still an exact sample
(`TK2.sample_path`), just a chosen one.
"""
function tk2_schematic_path(leading; nres=4, tau_min=1.2, tau_max=6.0)
    for seed in 1:10_000
        segments = TK2Model.sample_path(leading, 1.0, nres; rng=Xoshiro(seed))
        all(seg -> tau_min <= seg.t[end] - seg.t[1] <= tau_max, segments) && return segments
    end
    return error("no schematic PDMP path found")
end

# --- Correction figure: panels (a) to (d) -------------------------------------------

"""
The ε of the exact samples in panel (a): a regular logarithmic lattice anchored at the first
row ε. Ten intervals per decade therefore put ε = 0.02 and ε = 0.2 exactly ten intervals apart,
without the uneven gaps created by inserting them into an unrelated fixed-length grid.
"""
function tk2_exact_epsilons(settings::TK2Settings)
    (; exact_points_per_decade, epsilon_plot_min, epsilon_plot_max, epsilon_rows) = settings
    anchor = first(epsilon_rows)
    first_step = round(Int, exact_points_per_decade * log10(epsilon_plot_min / anchor))
    last_step = floor(Int, exact_points_per_decade * log10(epsilon_plot_max / anchor))
    epsilons = anchor .* 10.0 .^ ((first_step:last_step) ./ exact_points_per_decade)
    for epsilon in epsilon_rows
        index = argmin(abs.(log10.(epsilons ./ epsilon)))
        require(
            isapprox(epsilons[index], epsilon; rtol=1e-12),
            "row epsilon=$epsilon does not lie on the exact-sample lattice",
        )
        epsilons[index] = epsilon
    end
    require(
        all(isapprox.(diff(log10.(epsilons)), 1 / exact_points_per_decade; atol=1e-12)),
        "exact epsilon samples are not logarithmically regular",
    )
    return epsilons
end

"""Panel (a): the mean branch coordinate of the O(ε) PDMP and of the exact chain."""
function tk2_mean_curves(model, settings, exact)
    (; epsilon_plot_min, epsilon_plot_max, curve_points) = settings
    epsilon_theory =
        10 .^ range(log10(epsilon_plot_min), log10(epsilon_plot_max); length=curve_points)
    epsilon_exact = tk2_exact_epsilons(settings)
    return MeanCurves(;
        epsilon_theory,
        corrected_y=[
            TK2Model.mean_branch_coordinate(tk2_corrected(model, epsilon)) for
            epsilon in epsilon_theory
        ],
        epsilon_exact,
        exact_y=[exact(epsilon).mean_y for epsilon in epsilon_exact],
    )
end

"""
Panel (b): at each row ε, the smoothed PDMP densities on a grid around the leading support,
the exact lattice density, and the errors between them.
"""
function tk2_density_rows(model, settings, leading, exact)
    (; population, points) = settings
    cutoff = tk2_cutoff(model, settings)
    lower, upper = TK2Model.support(leading)
    grid = collect(
        range(
            max(lower - 2.5 / sqrt(population), 1e-3), upper + 3 / sqrt(population); length=points
        ),
    )
    smoothed(solution) = TK2Model.finite_size_profile(solution, population; grid).total
    relative_sup(theory, exact) = maximum(abs.(theory .- exact)) / maximum(exact)
    return map(settings.epsilon_rows) do epsilon
        corrected = tk2_corrected(model, epsilon)
        e = exact(epsilon)
        exact_density = e.masses .* population
        nlo = max(0, floor(Int, first(grid) * population))
        nhi = min(cutoff, ceil(Int, last(grid) * population))
        exact_on_grid = exact_density[clamp.(round.(Int, grid .* population) .+ 1, 1, cutoff + 1)]
        density = Comparison(exact_on_grid, smoothed(leading), smoothed(corrected))
        return DensityRow(;
            epsilon,
            grid,
            density,
            bar_edges=((nlo:(nhi + 1)) .- 0.5) ./ population,
            bar_heights=exact_density[(nlo + 1):(nhi + 1)],
            mean_y_exact=e.mean_y,
            mean_y_corrected=TK2Model.mean_branch_coordinate(corrected),
            l1=ByOrder(
                TK2Model.density_l1(leading, e.masses, population),
                TK2Model.density_l1(corrected, e.masses, population),
            ),
            sup=errors(density, relative_sup),
        )
    end
end

"""
Wasserstein-1 and lattice L1 errors of both smoothed PDMPs against the exact chain at each
`(ε, V)`, as `(w1, l1)`.
"""
function tk2_errors(model, leading, exact, points)
    function both(metric)
        error(solution, epsilon, V) = metric(solution, exact(epsilon, V).masses, V)
        return ByOrder(
            [error(leading, epsilon, V) for (epsilon, V) in points],
            [error(tk2_corrected(model, epsilon), epsilon, V) for (epsilon, V) in points],
        )
    end
    return both(TK2Model.density_w1), both(TK2Model.density_l1)
end

"""Panel (c): the L1 errors against ε at each V of `error_populations_epsilon`."""
function tk2_epsilon_sweeps(model, settings, leading, exact)
    return [
        EpsilonSweep(
            V,
            settings.error_epsilons,
            tk2_errors(
                model, leading, exact, [(epsilon, V) for epsilon in settings.error_epsilons]
            )...,
        ) for V in settings.error_populations_epsilon
    ]
end

"""
Panel (d): the L1 errors against μ = 1/V at each ε of `error_mu_epsilons`, small enough that
the ε error does not mask the finite-volume one.
"""
function tk2_mu_sweeps(model, settings, leading, exact)
    return [
        MuSweep(
            epsilon,
            1 ./ settings.error_populations_mu,
            tk2_errors(
                model, leading, exact, [(epsilon, V) for V in settings.error_populations_mu]
            )...,
        ) for epsilon in settings.error_mu_epsilons
    ]
end

# --- Pipeline -----------------------------------------------------------------------

function compute_tk2(model::TK2ModelConfig, settings::TK2Settings)
    leading = tk2_leading(model)
    exact, solved = tk2_exact_solver(model, settings)
    means = tk2_mean_curves(model, settings, exact)
    rows = tk2_density_rows(model, settings, leading, exact)
    epsilon_sweeps = tk2_epsilon_sweeps(model, settings, leading, exact)
    mu_sweeps = tk2_mu_sweeps(model, settings, leading, exact)
    return TK2Results(;
        leading,
        leading_y=TK2Model.mean_branch_coordinate(leading),
        cutoff=tk2_cutoff(model, settings),
        maximum_boundary_mass=maximum(e.boundary_mass for e in values(solved)),
        path=tk2_schematic_path(leading),
        means,
        rows,
        epsilon_sweeps,
        mu_sweeps,
    )
end

"""Fail unless the cutoff is safe and the errors show the expected O(ε) and O(ε²) regimes."""
function check_tk2(r::TK2Results)
    r.maximum_boundary_mass < 1e-8 || error("TK2 finite-state projection cutoff mass exceeds 1e-8")
    finest = last(r.epsilon_sweeps)
    check_improves(finest.w1, "the TK2 density error at the largest ε")
    maximum(local_slopes(finest.epsilons, finest.w1.leading)) > 0.7 ||
        error("TK2 leading error lacks the expected O(epsilon) regime")
    maximum(local_slopes(finest.epsilons, finest.w1.corrected)) > 1.3 ||
        error("TK2 corrected error lacks the expected O(epsilon^2) regime")
    return nothing
end

"""Data for the manuscript's `figures/tk2-leading.tex`: the schematic path and its residences."""
function write_tk2_leading(output, model, r::TK2Results)
    data = joinpath(output, "data", "tk2-leading")
    path = r.path
    fixed_points = TK2Model.support(r.leading)
    write_columns(
        joinpath(data, "path.csv"),
        ("t", "s", "plus"),
        reduce(vcat, [segment.t for segment in path]),
        reduce(vcat, [segment.x for segment in path]),
        reduce(vcat, [fill(segment.branch === :plus, length(segment.t)) for segment in path]),
    )
    write_columns(
        joinpath(data, "residences.csv"),
        ("t0", "t1", "plus", "s0", "s1"),
        [first(segment.t) for segment in path],
        [last(segment.t) for segment in path],
        [segment.branch === :plus for segment in path],
        [first(segment.x) for segment in path],
        [last(segment.x) for segment in path],
    )
    return write_tex_macros(
        joinpath(data, "params.tex"),
        (
            "fixedPointPlus" => fixed_points[1],
            "fixedPointMinus" => fixed_points[2],
            "pathEndTime" => last(last(path).t),
            "kappaValue" => model.kappa,
            "betaValue" => model.beta,
        ),
    )
end

"""
Data for the manuscript's `figures/tk2-correction.tex`: one CSV per curve, and the macros of
`tk2_correction_macros`.
"""
function write_tk2_correction(output, model, settings, r::TK2Results)
    data = joinpath(output, "data", "tk2-correction")
    for (row, density_row) in zip(("a", "b"), r.rows)
        (; grid, density, bar_edges, bar_heights) = density_row
        write_columns(
            joinpath(data, "density-$row.csv"),
            ("s", "leading", "corrected"),
            grid,
            density.leading,
            density.corrected,
        )
        # Left edges with the last bar repeated, for a const plot of the exact lattice density.
        write_columns(
            joinpath(data, "exact-$row.csv"),
            ("s", "density"),
            bar_edges,
            [bar_heights; bar_heights[end]],
        )
    end
    write_columns(
        joinpath(data, "mean-corrected.csv"),
        ("epsilon", "y"),
        r.means.epsilon_theory,
        r.means.corrected_y,
    )
    write_columns(
        joinpath(data, "mean-exact.csv"), ("epsilon", "y"), r.means.epsilon_exact, r.means.exact_y
    )
    one, two = r.epsilon_sweeps
    write_columns(
        joinpath(data, "errors-epsilon.csv"),
        ("epsilon", "leadingone", "correctedone", "leadingtwo", "correctedtwo"),
        one.epsilons,
        one.w1.leading,
        one.w1.corrected,
        two.w1.leading,
        two.w1.corrected,
    )
    one, two = r.mu_sweeps
    write_columns(
        joinpath(data, "errors-mu.csv"),
        ("mu", "leadingone", "correctedone", "leadingtwo", "correctedtwo"),
        one.mus,
        one.w1.leading,
        one.w1.corrected,
        two.w1.leading,
        two.w1.corrected,
    )
    return write_tex_macros(joinpath(data, "params.tex"), tk2_correction_macros(model, settings, r))
end

"""
Numbers that `figures/tk2-correction.tex` reads as TeX macros: the parameters it prints, axis
limits and guide anchors. The padding factors were tuned by eye for that figure. Change them
here, not in the .tex.
"""
function tk2_correction_macros(model, settings, r::TK2Results)
    (; rows, epsilon_sweeps, mu_sweeps) = r
    finest = last(epsilon_sweeps)
    mus = first(mu_sweeps).mus
    all_errors(sweeps) = reduce(vcat, [[s.w1.leading; s.w1.corrected] for s in sweeps])
    guide_gap = 2  # reference guides sit this factor below the curve they follow

    # (a) Mean y against ε. Leave 8% of the data span below the lowest mean, and 15% above the
    # leading y, for the legend. The ε axis, shared with (c), runs a factor 1.25 past the
    # plotted range.
    lowest_mean = min(minimum(r.means.exact_y), minimum(r.means.corrected_y))
    meanymin, meanymax = pad_linear(lowest_mean, r.leading_y; below=0.08, above=0.15)
    epsmin, epsmax = pad_log(settings.epsilon_plot_min, settings.epsilon_plot_max; below=1.25)
    # (b) Densities: 20% headroom over each row's exact peak, on the theory grid.
    density_headroom = 1.2
    # (c) Errors against ε: a factor 1.6 below the smallest error and 0.2 decade above the
    # largest. The O(ε) and O(ε²) guides start at the largest ε, below the finest-V curves.
    errcmin, errcmax = pad_log(extrema(all_errors(epsilon_sweeps))...; below=1.6, above=10^0.2)
    # (d) Errors against μ: from the decade below the smallest error to 40% above the largest.
    # The guide has the fitted finite-volume slope of the corrected errors (mean over both ε)
    # and starts at the largest μ, below the last corrected curve.
    mu_errors = all_errors(mu_sweeps)
    mumin, mumax = pad_log(extrema(mus)...; below=1.25)

    return (
        "kappaValue" => model.kappa,
        "betaValue" => model.beta,
        "densityVolume" => settings.population,
        "leadingMeanY" => r.leading_y,
        "meanYMin" => meanymin,
        "meanYMax" => meanymax,
        "epsAxisMin" => epsmin,
        "epsAxisMax" => epsmax,
        "epsRowA" => rows[1].epsilon,
        "epsRowB" => rows[2].epsilon,
        "densityMaxA" => density_headroom * maximum(rows[1].density.exact),
        "densityMaxB" => density_headroom * maximum(rows[2].density.exact),
        "sAxisMin" => first(rows[1].grid),
        "sAxisMax" => last(rows[1].grid),
        "errorVolumeOne" => epsilon_sweeps[1].population,
        "errorVolumeTwo" => epsilon_sweeps[2].population,
        "epsErrorMin" => errcmin,
        "epsErrorMax" => errcmax,
        "muErrorEpsOne" => mu_sweeps[1].epsilon,
        "muErrorEpsTwo" => mu_sweeps[2].epsilon,
        "muAxisMin" => mumin,
        "muAxisMax" => mumax,
        "muErrorMin" => decade_below(minimum(mu_errors)),
        "muErrorMax" => 1.4 * maximum(mu_errors),
        "epsGuideX" => last(finest.epsilons),
        "leadingGuideY" => last(finest.w1.leading) / guide_gap,
        "correctedGuideY" => last(finest.w1.corrected) / guide_gap,
        "muGuideSlope" =>
            sum(loglog_slope(mus, s.w1.corrected) for s in mu_sweeps) / length(mu_sweeps),
        "muGuideX" => maximum(mus),
        "muGuideY" => last(mu_sweeps).w1.corrected[argmax(mus)] / guide_gap,
    )
end

"""
Numbers for the captions and text, in tables that mirror `TK2Results`: the panel (b) rows, the
sweeps of (c, d) with their log-log slopes, and the cutoff mass.
"""
function tk2_metrics(r::TK2Results)
    (; epsilon_sweeps, mu_sweeps) = r
    function row_table(row)
        return Dict{String,Any}(
            "mean_y_exact" => row.mean_y_exact,
            "mean_y_corrected" => row.mean_y_corrected,
            "l1" => metrics_table(row.l1),
            "sup" => metrics_table(row.sup),
        )
    end
    function epsilon_table(s)
        return Dict{String,Any}(
            "w1" => metrics_table(s.w1),
            "l1" => metrics_table(s.l1),
            "local_slopes" => metrics_table(map(e -> local_slopes(s.epsilons, e), s.w1)),
        )
    end
    function mu_table(s)
        return Dict{String,Any}(
            "w1" => metrics_table(s.w1),
            "l1" => metrics_table(s.l1),
            "slope" => metrics_table(map(e -> loglog_slope(s.mus, e), s.w1)),
        )
    end
    return Dict{String,Any}(
        "cutoff" => r.cutoff,
        "leading_mean_y" => r.leading_y,
        "maximum_boundary_mass" => r.maximum_boundary_mass,
        "schematic_switch_states" => [last(segment.x) for segment in r.path[1:(end - 1)]],
        "density_rows" => Dict(epsilon_key(row.epsilon) => row_table(row) for row in r.rows),
        "error_epsilon_sweeps" =>
            Dict("V$(s.population)" => epsilon_table(s) for s in epsilon_sweeps),
        "error_mu_sweeps" => Dict(epsilon_key(s.epsilon) => mu_table(s) for s in mu_sweeps),
        # At small ε the corrected error is the finite-volume error alone, the same at each ε.
        "error_mu_sweep_corrected_relative_spread" =>
            maximum(abs.(first(mu_sweeps).w1.corrected ./ last(mu_sweeps).w1.corrected .- 1)),
    )
end

"""
Reproduce the data of the manuscript's two TK2 figures in `output/<mode>/tk2/`: compute
everything, fail unless the results show what the figures claim, then write the figure data
and `metrics.toml`.
"""
function reproduce_tk2(mode::Symbol)
    config = load_config(TK2Config, "tk2")
    model, settings = config.model, getfield(config, mode)
    output = output_dir("tk2", mode)
    results = compute_tk2(model, settings)
    check_tk2(results)
    write_tk2_leading(output, model, results)
    write_tk2_correction(output, model, settings, results)
    save_metrics(output, mode, model, settings, tk2_metrics(results))
    return output
end
