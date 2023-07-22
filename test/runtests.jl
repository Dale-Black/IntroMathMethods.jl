using IntroMathMethods
using Test
using Symbolics, SymbolicNumericIntegration

@testset "IntroMathMethods.jl" begin
    @testset "Integrals" begin
        @variables x
        f = exp(-x)
        @test integrate_definite(f, x, 0, Inf) == 1
    end

    @testset "Taylor Series" begin
        @variables x
        f = sin(x) / (cosh(x) + 2)
        f_ts = taylor_series(f, x, 3)
        @test isapprox(substitute(f, x => 0.1), substitute(f_ts, x => 0.1); atol = 1e-3)
        @test isapprox(substitute(f, x => -0.1), substitute(f_ts, x => -0.1); atol = 1e-3)

        exp_ts = taylor_series_exp(x, 10)
        @test isapprox(exp(1), Float64(substitute(exp_ts, x => 1).val); atol = 1e-2)
        @test isapprox(exp(2), Float64(substitute(exp_ts, x => 2).val); atol = 1e-2)
    end

    @testset "Lagrange Multipliers" begin
        @variables x y
        f = x^2 + y^2
        g = x + 2y - 4	
        x0, y0, _ = lagrangian(f, g, x, y)
        @test (x0, y0) == (4/5, 8/5)
    end
end
