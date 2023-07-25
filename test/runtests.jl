using IntroMathMethods
using PythonCall
using Test

@testset "IntroMathMethods.jl" begin
    @testset "Lagrangian" begin
        x, y, λ = sp.symbols("x, y, λ")
        f = x^2 + y^2
        g = x + 2y - 4
        l = lagrangian(f, g, x, y)
        dict = pyconvert(Dict, l)
        @test (haskey(dict, x) && haskey(dict, y) && haskey(dict, λ))
        @test (8//5 in values(dict) && 4//5 in values(dict))
    end
end
