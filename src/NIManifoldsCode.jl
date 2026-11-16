module NIManifoldsCode

using SparseArrays

"""
    stationary_vector(generator) -> Vector{Float64}

Stationary probability vector π of a continuous-time Markov chain whose sparse `generator` has
rows = from-state and zero row sums, so that π Q = 0. The first balance equation is replaced
by the normalisation Σπ = 1 and the rest are solved by sparse LU. Shared by the exact chains of
`TK2` and `Ants`.
"""
function stationary_vector(generator::SparseMatrixCSC)
    A = sparse(transpose(generator))
    x = A[2:end, 2:end] \ Vector(-A[2:end, 1])
    stationary = vcat(1.0, x)
    return stationary ./ sum(stationary)
end

include("tk2/TK2.jl")
include("ants/Ants.jl")

export TK2, Ants

end
