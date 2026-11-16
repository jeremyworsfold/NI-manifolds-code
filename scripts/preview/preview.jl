#!/usr/bin/env julia
#
# Makie quick-looks of the manuscript figures, drawn from the figure data that
# `scripts/reproduce.jl` writes: the same CSV files and `params.tex` macros the manuscript's
# pgfplots figures read. Optional, and in its own environment so the reproduction never loads
# Makie:
#
#     julia --project=scripts/preview scripts/preview/preview.jl --mode quick --example all
#
# Each preview PDF is written next to its example's `metrics.toml`.

using DelimitedFiles

include(joinpath(@__DIR__, "plotting.jl"))
set_theme!(include(joinpath(@__DIR__, "theme.jl")))
CairoMakie.activate!()

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const EXAMPLES = ("tk2", "ants")

"""Columns of a figure-data CSV as a NamedTuple keyed by its header (`nan` reads as NaN)."""
function read_columns(path)
    data, header = readdlm(path, ','; header=true)
    return NamedTuple(Symbol(h) => Float64.(data[:, i]) for (i, h) in enumerate(header))
end

"""The `\\def\\name{value}` macros of a figure's `params.tex`, as numbers keyed by name."""
function read_macros(path)
    text = read(path, String)
    return Dict(m[1] => parse(Float64, m[2]) for m in eachmatch(r"\\def\\(\w+)\{([^}]*)\}", text))
end

"""A figure's data directory, its columns reader and its macros."""
function figure_data(output, figure)
    data = joinpath(output, "data", figure)
    return (name -> read_columns(joinpath(data, "$name.csv"))),
    read_macros(joinpath(data, "params.tex"))
end

panel_label!(position, label) = Label(position, label; font=:bold, padding=(0, 6, 2, 0))

const EXACT_MARKER = (; color=EXACT_COLOR, strokewidth=0.5, strokecolor=:black)

"""Legend entries for the three roles every figure shares: O(1), O(ε) and the exact chain."""
function role_legend()
    return (
        [
            LineElement(; color=LEADING_COLOR, linewidth=1.0),
            LineElement(; color=CORRECTED_COLOR, linewidth=1.8),
            [
                PolyElement(; color=(EXACT_COLOR, 0.5)),
                MarkerElement(; EXACT_MARKER..., marker=:circle, markersize=6),
            ],
        ],
        [L"\mathcal{O}(1)", L"\mathcal{O}(\varepsilon)", L"\mathrm{Exact}"],
    )
end

# --- TK2 -------------------------------------------------------------------------------

const TK2_BRANCH_COLORS = Dict(true => COLORS[4], false => COLORS[5])  # plus, minus

"""The schematic path as one `(t, s, plus)` segment per residence."""
function tk2_segments(path)
    breaks = [0; findall(diff(path.plus) .!= 0); length(path.t)]
    return [
        (; t=path.t[(a + 1):b], s=path.s[(a + 1):b], plus=path.plus[b] == 1) for
        (a, b) in zip(breaks[1:(end - 1)], breaks[2:end])
    ]
end

"""Leading-order PDMP in the (x₊, x₋) plane: flow along the manifolds, jumps along fibres."""
function draw_tk2_plane!(position, segments, fixed_points)
    limit, inset = 1.95, 0.08  # the inset keeps the thick manifold strokes off the spines
    axis = Axis(
        position;
        xlabel=L"x_+",
        ylabel=L"x_-",
        aspect=DataAspect(),
        limits=(-inset, limit, -inset, limit),
        xticks=[0, 1],
        yticks=[0, 1],
    )
    for s in 0.25:0.25:(2limit)
        lines!(axis, [0, s], [s, 0]; color=(:black, 0.15), linewidth=0.4, linestyle=:dash)
    end
    for segment in segments
        xy = segment.plus ? (segment.s, zero(segment.s)) : (zero(segment.s), segment.s)
        lines!(axis, xy...; color=TK2_BRANCH_COLORS[segment.plus], linewidth=2.5)
    end
    for segment in segments[1:(end - 1)]
        s = last(segment.s)
        lines!(axis, [s, 0], [0, s]; color=:black, linewidth=0.6, linestyle=:dot)
    end
    plus, minus = fixed_points
    scatter!(
        axis, [plus, 0], [0, minus]; color=:white, strokecolor=:black, strokewidth=0.8, markersize=6
    )
    text!(axis, plus + 0.08, 0.05; text=L"s^*_+", align=(:left, :bottom))
    text!(axis, 0.06, minus; text=L"s^*_-", align=(:left, :center))
    return axis
end

"""The same path as s(t), coloured by branch."""
function draw_tk2_trace!(position, segments, fixed_points)
    axis = Axis(position; xlabel=L"t", ylabel=L"s")
    hlines!(
        axis,
        collect(fixed_points);
        color=[TK2_BRANCH_COLORS[true], TK2_BRANCH_COLORS[false]],
        linestyle=:dash,
        linewidth=0.6,
    )
    for segment in segments
        lines!(axis, segment.t, segment.s; color=TK2_BRANCH_COLORS[segment.plus], linewidth=1.4)
    end
    vlines!(
        axis,
        [last(segment.t) for segment in segments[1:(end - 1)]];
        color=(:black, 0.3),
        linewidth=0.5,
        linestyle=:dot,
    )
    xlims!(axis, 0, last(last(segments).t))
    return axis
end

"""Ticks at 1, 2 and 5 times each decade between `lo` and `hi`, labelled as decimals."""
function ticks_125(lo, hi)
    ticks = [m * 10.0^d for d in floor(Int, log10(lo)):ceil(Int, log10(hi)) for m in (1, 2, 5)]
    ticks = filter(t -> lo <= t <= hi, ticks)
    return (ticks, [string(round(t; sigdigits=1)) for t in ticks])
end

"""
Log-log density error against a swept `x`, from an `errors-*.csv` table. Colour marks the order
and line style the two sweeps, dashed (`one`) then solid (`two`). Each guide
`(power, label, x0, y0)` is a faint line `y0 (x/x0)^power` over the decade below `x0`.
"""
function draw_tk2_error_panel!(
    position, x, errors, labels; xlabel, xlimits, ylimits, guides, legend_position
)
    decades = ceil(Int, log10(ylimits[1])):floor(Int, log10(ylimits[2]))
    axis = Axis(
        position;
        xlabel,
        ylabel=L"\ell^1\ \mathrm{error}",
        xscale=log10,
        yscale=log10,
        xticks=ticks_125(xlimits...),
        yticks=(10.0 .^ decades, [L"10^{%$d}" for d in decades]),
        yminorticks=IntervalsBetween(9),
        yminorticksvisible=true,
    )
    for (power, label, x0, y0) in guides
        span = [x0 / 10, x0]
        lines!(axis, span, y0 .* (span ./ x0) .^ power; color=(:black, 0.4), linewidth=0.7)
        middle = x0 / sqrt(10)
        text!(
            axis,
            middle,
            y0 * (middle / x0)^power;
            text=label,
            fontsize=8,
            color=(:black, 0.6),
            align=(:left, :top),
            offset=(2, -2),
        )
    end
    linestyles = (:dash, :solid)
    for (sweep, linestyle) in zip(("one", "two"), linestyles)
        for (order, color) in (("leading", LEADING_COLOR), ("corrected", CORRECTED_COLOR))
            values = errors[Symbol(order, sweep)]
            scatterlines!(axis, x, values; color, linestyle, linewidth=1.3, markersize=4)
        end
    end
    limits!(axis, xlimits, ylimits)
    axislegend(
        axis,
        [LineElement(; color=(:black, 0.7), linestyle, linewidth=1.3) for linestyle in linestyles],
        labels;
        position=legend_position,
        orientation=:horizontal,
        padding=3,
        colgap=8,
        patchsize=(18, 7),
        patchlabelgap=4,
        framevisible=false,
        labelsize=8,
    )
    return axis
end

function preview_tk2(output)
    columns, macros = figure_data(output, "tk2-leading")
    segments = tk2_segments(columns("path"))
    fixed_points = (macros["fixedPointPlus"], macros["fixedPointMinus"])
    schematic = fig_in_cm(5.8, 17.0)
    draw_tk2_plane!(schematic[1, 1], segments, fixed_points)
    draw_tk2_trace!(schematic[1, 2], segments, fixed_points)
    colsize!(schematic.layout, 1, Aspect(1, 1.0))
    colgap!(schematic.layout, Relative(0.06))
    panel_label!(schematic[1, 1, TopLeft()], "a")
    panel_label!(schematic[1, 2, TopLeft()], "b")
    save_fig("TK-leading-order-preview", schematic; path=output)

    columns, m = figure_data(output, "tk2-correction")
    figure = fig_in_cm(12.6, 17.0)
    left = figure[1, 1] = GridLayout()
    right = figure[1, 2] = GridLayout()
    epsilon_limits = (m["epsAxisMin"], m["epsAxisMax"])

    # (a) Mean branch coordinate against ε, with guides at the ε of the density rows.
    mean_axis = Axis(
        left[1, 1]; ylabel=L"\langle y\rangle", xscale=log10, xticks=ticks_125(epsilon_limits...)
    )
    vlines!(mean_axis, [m["epsRowA"], m["epsRowB"]]; color=(:black, 0.22), linewidth=0.5)
    corrected, exact = columns("mean-corrected"), columns("mean-exact")
    lines!(mean_axis, corrected.epsilon, corrected.y; color=CORRECTED_COLOR, linewidth=1.8)
    hlines!(mean_axis, m["leadingMeanY"]; color=LEADING_COLOR, linewidth=1.0)
    scatter!(mean_axis, exact.epsilon, exact.y; EXACT_MARKER..., markersize=6)
    limits!(mean_axis, epsilon_limits, (m["meanYMin"], m["meanYMax"]))
    axislegend(
        mean_axis,
        role_legend()...;
        position=:lb,
        padding=3,
        rowgap=1,
        patchsize=(18, 8),
        patchlabelgap=5,
        framevisible=false,
    )

    # (b) Stationary densities of s at the two row ε.
    for (row, (name, epsilon)) in enumerate((("a", m["epsRowA"]), ("b", m["epsRowB"])))
        density, exact = columns("density-$name"), columns("exact-$name")
        last_row = row == 2
        axis = Axis(right[row, 1]; ylabel=L"P_0^*(s)", xlabel=last_row ? L"s" : "")
        manual_hist!(axis, exact.s, exact.density[1:(end - 1)]; color=EXACT_COLOR)
        lines!(axis, density.s, density.corrected; color=CORRECTED_COLOR, linewidth=1.8)
        lines!(axis, density.s, density.leading; color=LEADING_COLOR, linewidth=1.0)
        text!(axis, 0.03, 0.86; text=L"\varepsilon = %$(epsilon)", space=:relative)
        limits!(axis, (m["sAxisMin"], m["sAxisMax"]), (0, m["densityMax" * uppercase(name)]))
        last_row || hidexdecorations!(axis; grid=false)
    end
    rowgap!(right, 6)

    # (c) Errors against ε at two V, with O(ε) and O(ε²) guides.
    errors = columns("errors-epsilon")
    V1, V2 = Int(m["errorVolumeOne"]), Int(m["errorVolumeTwo"])
    epsilon_axis = draw_tk2_error_panel!(
        left[2, 1],
        errors.epsilon,
        errors,
        [L"V = %$V1", L"V = %$V2"];
        xlabel=L"\varepsilon",
        xlimits=epsilon_limits,
        ylimits=(m["epsErrorMin"], m["epsErrorMax"]),
        guides=(
            (1, L"\varepsilon", m["epsGuideX"], m["leadingGuideY"]),
            (2, L"\varepsilon^2", m["epsGuideX"], m["correctedGuideY"]),
        ),
        legend_position=(0.03, 0.99),
    )

    # (d) Errors against μ = 1/V at two small ε, with the fitted finite-volume slope.
    errors = columns("errors-mu")
    slope = round(m["muGuideSlope"]; digits=2)
    ε1, ε2 = m["muErrorEpsOne"], m["muErrorEpsTwo"]
    draw_tk2_error_panel!(
        left[3, 1],
        errors.mu,
        errors,
        [L"\varepsilon = %$ε1", L"\varepsilon = %$ε2"];
        xlabel=L"\mu",
        xlimits=(m["muAxisMin"], m["muAxisMax"]),
        ylimits=(m["muErrorMin"], m["muErrorMax"]),
        guides=((m["muGuideSlope"], L"\mu^{%$slope}", m["muGuideX"], m["muGuideY"]),),
        legend_position=(0.97, 0.03),
    )

    # (a) and (c) share the ε axis, so (a) drops its tick labels.
    linkxaxes!(mean_axis, epsilon_axis)
    hidexdecorations!(mean_axis; grid=false, ticks=false, minorticks=false)
    colsize!(figure.layout, 1, Relative(0.47))
    colgap!(figure.layout, Relative(0.05))
    rowsize!(left, 1, Relative(0.26))
    rowgap!(left, 8)
    panel_label!(left[1, 1, TopLeft()], "a")
    panel_label!(right[1, 1, TopLeft()], "b")
    panel_label!(left[2, 1, TopLeft()], "c")
    panel_label!(left[3, 1, TopLeft()], "d")
    return save_fig("TK-correction-preview", figure; path=output)
end

# --- Ant foraging ----------------------------------------------------------------------

"""The ant-foraging density and residence-time data as a compact preview."""
function preview_ants(output)
    name = "ants-correction"
    columns, m = figure_data(output, name)
    epsilon = m["epsPanel"]
    density, exact, survival = columns("density"), columns("exact"), columns("survival")
    figure = fig_in_cm(6.5, 17.0)

    # (a) Stationary food remaining at source A.
    density_axis = Axis(figure[1, 1]; xlabel=L"\mathrm{food\ level}\ u", ylabel=L"\tilde{\rho}(u)")
    edges = [exact.u; last(exact.u) + 1 / m["foodSteps"]]
    manual_hist!(density_axis, edges, exact.density; color=EXACT_COLOR)
    lines!(density_axis, density.u, density.corrected; color=CORRECTED_COLOR, linewidth=1.8)
    lines!(density_axis, density.u, density.leading; color=LEADING_COLOR, linewidth=1.0)
    limits!(density_axis, (-0.02, 1.08), (0, m["densityMax"]))
    axislegend(
        density_axis,
        role_legend()...;
        position=(0.45, 0.95),
        padding=3,
        rowgap=1,
        patchsize=(18, 8),
        patchlabelgap=5,
        framevisible=false,
    )
    text!(
        density_axis,
        0.45,
        0.62;
        text=L"\varepsilon = %$epsilon",
        space=:relative,
        align=(:center, :top),
    )

    # (b) Line style marks the source; exact values use circles and squares.
    sources = (
        (; key="A", name="source A", linestyle=:solid, marker=:circle),
        (; key="B", name="source B", linestyle=:dash, marker=:rect),
    )
    decades = -2:0
    survival_axis = Axis(
        figure[1, 2];
        yscale=log10,
        ylabel=L"S(t)",
        xlabel=L"t",
        yticks=(10.0 .^ decades, [L"10^{%$d}" for d in decades]),
    )
    times = survival.t
    pick = 1:max(1, length(times) ÷ 20):length(times)
    for (; key, linestyle, marker) in sources
        corrected = survival[Symbol("corrected", key)]
        leading = survival[Symbol("leading", key)]
        exact = survival[Symbol("exact", key)]
        lines!(survival_axis, times, corrected; color=CORRECTED_COLOR, linewidth=1.8, linestyle)
        lines!(survival_axis, times, leading; color=LEADING_COLOR, linewidth=1.0, linestyle)
        scatter!(survival_axis, times[pick], exact[pick]; EXACT_MARKER..., marker, markersize=5)
    end
    axislegend(
        survival_axis,
        [
            [
                LineElement(; color=(:black, 0.7), b.linestyle, linewidth=1.3),
                MarkerElement(; EXACT_MARKER..., b.marker, markersize=5),
            ] for b in sources
        ],
        [b.name for b in sources];
        position=:lb,
        padding=3,
        rowgap=1,
        patchsize=(18, 7),
        patchlabelgap=4,
        framevisible=false,
        labelsize=8,
    )
    limits!(survival_axis, (0, m["survivalSpan"]), (0.8 * 10.0^first(decades), 1.3))

    colgap!(figure.layout, Relative(0.04))
    panel_label!(figure[1, 1, TopLeft()], "a")
    panel_label!(figure[1, 2, TopLeft()], "b")
    return save_fig("$name-preview", figure; path=output)
end

# --- Command line ----------------------------------------------------------------------

"""The value after `flag` in `arguments`, or `default` when the flag is absent."""
function option(arguments, flag, default)
    i = findfirst(==(flag), arguments)
    isnothing(i) && return default
    i < length(arguments) || error("$flag needs a value")
    return arguments[i + 1]
end

function main(arguments=ARGS)
    mode = option(arguments, "--mode", "quick")
    example = option(arguments, "--example", "all")
    mode in ("quick", "full") || error("--mode must be quick or full")
    example in ("all", EXAMPLES...) ||
        error("--example must be one of: all, $(join(EXAMPLES, ", "))")
    for item in (example == "all" ? EXAMPLES : (example,))
        output = joinpath(REPOSITORY_ROOT, "output", mode, item)
        isdir(joinpath(output, "data")) ||
            error("no figure data in $output: run scripts/reproduce.jl --mode $mode first")
        item == "tk2" ? preview_tk2(output) : preview_ants(output)
        println("drew $item previews in $output")
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
