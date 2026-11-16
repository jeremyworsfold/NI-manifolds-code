# The vocabulary both examples compare in: the exact chain, the leading-order PDMP and the O(ε)
# PDMP, and how errors between them converge.

"""
A quantity from the leading-order PDMP and from the O(ε) PDMP, most often their errors against
the exact chain.
"""
struct ByOrder{T}
    leading::T
    corrected::T
end

"""The same quantity from the exact chain, the leading-order PDMP and the O(ε) PDMP."""
struct Comparison{T}
    exact::T
    leading::T
    corrected::T
end

"""Distances `distance(theory, exact)` of both PDMPs from the exact chain."""
function errors(c::Comparison, distance)
    ByOrder(distance(c.leading, c.exact), distance(c.corrected, c.exact))
end

"""`f` applied to both orders."""
Base.map(f, x::ByOrder) = ByOrder(f(x.leading), f(x.corrected))

"""
Throw unless the O(ε) error is below the leading-order one. For errors along a sweep, compare the
last point, the largest ε or μ, where the correction matters most.
"""
function check_improves(errors::ByOrder, what::AbstractString)
    final(x) = x isa AbstractVector ? last(x) : x
    final(errors.corrected) < final(errors.leading) ||
        error("the O(ε) correction does not improve $what")
    return nothing
end

"""A `[leading, corrected]` table for `metrics.toml`."""
metrics_table(x::ByOrder) = Dict{String,Any}("leading" => x.leading, "corrected" => x.corrected)

"""Least-squares slope of `log y` against `log x`."""
function loglog_slope(x, y)
    lx, ly = log.(x), log.(y)
    dx = lx .- sum(lx) / length(lx)
    return sum(dx .* (ly .- sum(ly) / length(ly))) / sum(abs2, dx)
end

"""Log-log slopes between neighbouring points, to read off O(ε), O(ε²) and the floor."""
local_slopes(x, y) = diff(log.(y)) ./ diff(log.(x))
