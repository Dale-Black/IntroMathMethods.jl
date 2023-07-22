module IntroMathMethods

using Symbolics, SymbolicNumericIntegration
using Symbolics: derivative, solve_for

include("series.jl")
include("lagrange.jl")
include("integrals.jl")

end