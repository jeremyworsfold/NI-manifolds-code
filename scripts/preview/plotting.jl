using CairoMakie
using LaTeXStrings

# ColorBrewer Dark2, in its published order.
const DARK2_COLORS = [
    "#1B9E77", # teal
    "#D95F02", # orange
    "#7570B3", # purple
    "#E7298A", # magenta
    "#66A61E", # green
    "#E6AB02", # mustard
    "#A6761D", # brown
    "#666666", # grey
]
const COLORS = DARK2_COLORS

# The first three colours are roles shared by every figure. Other categories, such as the two
# branches or two parameter values, take colours from the rest of the palette.
const EXACT_COLOR = COLORS[1]
const LEADING_COLOR = COLORS[2]
const CORRECTED_COLOR = COLORS[3]

const PIXEL_SCALE = 28.3465

"""Create a figure whose `height` and `width` are specified in centimetres."""
function fig_in_cm(height::Real, width::Real; padding=(2.5, 5, 2.5, 5), args...)
    return Figure(;
        size=(width * PIXEL_SCALE, height * PIXEL_SCALE), figure_padding=padding, args...
    )
end

function save_fig(name::AbstractString, figure::Figure; path="figures", fmt="pdf", args...)
    mkpath(path)
    destination = joinpath(path, "$name.$fmt")
    CairoMakie.save(destination, figure; pt_per_unit=1, args...)
    return destination
end

"""
Fill a lattice density as one step region without internal bin lines: bar `i` spans
`edges[i]` to `edges[i + 1]` at height `heights[i]`.
"""
function manual_hist!(axis::Axis, edges, heights; color=COLORS[1], args...)
    xs = repeat(edges; inner=2)[2:(end - 1)]
    ys = repeat(heights; inner=2)
    return band!(axis, xs, zeros(length(ys)), ys; color=(color, 0.5), args...)
end
