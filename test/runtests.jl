using IntroMathMethods
using Test
using Symbolics

@testset "IntroMathMethods.jl" begin
    @variables x
    @testset "Taylor Series" begin
        f = sin(x) / (cosh(x) + 2)
        f_ts = taylor_series(f, x, 3)
        @test isapprox(substitute(f, x => 0.1), substitute(f_ts, x => 0.1); atol = 1e-3)
        @test isapprox(substitute(f, x => -0.1), substitute(f_ts, x => -0.1); atol = 1e-3)

        exp_ts = taylor_series_exp(x, 10)
        @test isapprox(exp(1), Float64(substitute(exp_ts, x => 1).val); atol = 1e-2)
        @test isapprox(exp(2), Float64(substitute(exp_ts, x => 2).val); atol = 1e-2)
    end
end
