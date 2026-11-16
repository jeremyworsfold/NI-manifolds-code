using NIManifoldsCode
using Integrals
using LinearAlgebra
using Random: Xoshiro
using SpecialFunctions: beta_inc, logbeta
using Test

const TK2Model = NIManifoldsCode.TK2
const AntsModel = NIManifoldsCode.Ants

include("tk2.jl")
include("ants.jl")
