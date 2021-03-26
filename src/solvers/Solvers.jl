module Solvers

using LinearAlgebra
using NLopt, NLsolve, NLSolvers,Roots
using  DiffResults, ForwardDiff
using StaticArrays
using PolynomialRoots
include("tunneling.jl")
include("ADNewton.jl")
include("fixpoint.jl")

polyroots(x) = PolynomialRoots.roots(x,polish=true)

end # module
