using IntroMathMethods
using Test
using Symbolics

@testset "IntroMathMethods.jl" begin
    @variables x
    @testset "Taylor Series" begin
        sym_exp = taylor_series_exp(x, 10)
        @test isapprox(exp(1), Float64(substitute(sym_exp, x => 1).val); atol = 1e-2)
        @test isapprox(exp(2), Float64(substitute(sym_exp, x => 2).val); atol = 1e-2)
    end
end
