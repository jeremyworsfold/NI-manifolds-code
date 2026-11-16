# Figure data for the manuscript's pgfplots figures: CSV columns and TeX macros.

"""How `write_columns` and `write_tex_macros` print a value: 8 significant digits, `nan` for NaN."""
format_value(x) = string(x)
format_value(x::Bool) = x ? "1" : "0"
format_value(x::Integer) = string(x)
format_value(x::Real) = isnan(x) ? "nan" : string(round(Float64(x); sigdigits=8))

"""
Write equal-length `columns` under `header` as comma-separated values, for the manuscript's
pgfplots figures. NaN is written as `nan`, which pgfplots reads as an unbounded point.
"""
function write_columns(path::AbstractString, header, columns...)
    length(header) == length(columns) || throw(ArgumentError("one header entry per column"))
    n = length(first(columns))
    all(length(column) == n for column in columns) ||
        throw(ArgumentError("columns differ in length"))
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        for i in 1:n
            println(io, join((format_value(column[i]) for column in columns), ","))
        end
    end
    return path
end

# --- Axis padding for the macros ----------------------------------------------------

"""Widen `[lo, hi]` by fractions of its span: `below` under `lo` and `above` over `hi`."""
pad_linear(lo, hi; below, above) = (lo - below * (hi - lo), hi + above * (hi - lo))

"""Widen a log-axis range `[lo, hi]` by the factor `below` under `lo` and `above` over `hi`."""
pad_log(lo, hi; below, above=below) = (lo / below, hi * above)

"""The power of ten at or below `x`, for a log axis that starts on a whole decade."""
decade_below(x) = 10.0^floor(log10(x))

"""
Write `\\def\\name{value}` lines for the numbers a pgfplots figure needs (limits, fixed points,
guide anchors). Macro names must be letters only.
"""
function write_tex_macros(path::AbstractString, macros)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "% Written by NI-manifolds-code (scripts/reproduce.jl); do not edit.")
        for (name, value) in macros
            all(isletter, name) || throw(ArgumentError("TeX macro names must be letters: $name"))
            println(io, "\\def\\$name{$(format_value(value))}")
        end
    end
    return path
end
