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
- `x`: Symbolic variable representing the x variable of `f`
- `y`: Symbolic variable representing the y variable of `f`

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
function lagrangian(f, g, x, y)
	@variables λ
	L = f - λ*g
	
	dLdx = derivative(L, x)
	dLdy = derivative(L, y)
	dLdλ = derivative(L, λ)

	sys = [dLdx ~ 0, dLdy ~ 0, dLdλ ~ 0]
	
	return solve_for(sys, [x, y, λ])
end

export lagrangian