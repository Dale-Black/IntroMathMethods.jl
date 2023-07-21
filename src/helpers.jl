using Symbolics, SymbolicNumericIntegration
using Symbolics: solve_for, derivative

export solve_for, derivative

"""
    taylor_series(f, x, n)

Compute the Taylor series expansion of a symbolic function `f` up to the `n`-th term, about the origin (0).

#### Arguments
- `f`: A symbolic expression representing a function. For example, if you have defined a symbolic variable `x` with `@syms x`, then `f` could be `sin(x)`.
- `x`: A symbolic variable representing the point about which the Taylor series is to be expanded. This variable is usually defined using the `@syms` macro.
- `n`: A non-negative integer representing the degree of the Taylor series expansion. The function will compute the Taylor series up to the `x^n` term.

#### Returns
- The function returns a symbolic expression representing the Taylor series expansion of `f` up to the `x^n` term.

#### Examples
```julia
@syms x
f = sin(x)
taylor_series(f, x, 5)
```
"""
function taylor_series(f, x, n::Int)
    f_ts = substitute(f, x => 0)
    for i in 1:n
		Dxi = Dx^i
        df = Dxi(f)
		df_0 = substitute(expand_derivatives(df), x => 0)
        f_ts += (df_0 * x^i) / factorial(i)
    end
    return f_ts
end

"""
    lagrangian(f, g)

Solve a constrained optimization problem using the Lagrangian method.

Given an objective function `f` and a constraint `g`, find the critical points 
that optimize `f` subject to `g` using the Lagrangian approach.

The Lagrangian `L` is formed as:

L(x, y, λ) = f(x, y) - λ(g(x, y))

The optimal points are found by taking derivatives of `L` w.r.t. x, y, λ 
and solving the resulting system of equations.

# Arguments
- `f`: Objective function to optimize
- `g`: Constraint equation 

# Returns
- Optimal values for x, y, λ that maximize/minimize `f` subject to `g`

# Examples
```julia
using Symbolics

@variables x y λ
f = x + y
g = x^2 + y^2 - 1
lagrangian(f, g)
```
"""
function lagrangian(f, g)
	@variables λ
	L = f - λ*g
	
	dLdx = derivative(L, x)
	dLdy = derivative(L, y)
	dLdλ = derivative(L, λ)

	sys = [dLdx ~ 0, dLdy ~ 0, dLdλ ~ 0]
	
	return solve_for(sys, [x, y, λ])
end

"""
    ∫(f, x)

Simple wrapper around `SymbolicNumericIntegration.integrate` to allow for the use of the integral symbol `∫` in Julia.
"""
function ∫(f, x)
	integrate(f, x)
end

function ∫(f)
	integrate(f)
end

"""
	integrate(f, x, x0, x1)

Simple wrapper around `integrate()` that allows for definite integrals to be computed.
"""
function SymbolicNumericIntegration.integrate(f, x, x0, x1)
	F = integrate(f, x)[1]
	F = substitute(F, x => x1) - substitute(F, x => x0)
end