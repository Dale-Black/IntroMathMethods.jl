"""
    lagrangian(f, g, x, y)

Solve a constrained optimization problem using the Lagrangian method.

Given an objective function `f` and a constraint `g`, find the critical points 
that optimize `f` subject to `g` using the Lagrangian approach.

The Lagrangian `L` is formed as:

``L(x, y, λ) = f(x, y) - λ*g(x, y)``

The optimal points are found by taking derivatives of `L` w.r.t. x, y, λ 
and solving the resulting system of equations.

# Arguments
- `f`: Objective function to optimize. This is a sympy expression.
- `g`: Constraint equation. This is a sympy expression.
- `x`: Symbolic variable representing the x variable of `f` and `g`.
- `y`: Symbolic variable representing the y variable of `f` and `g`.

# Returns
- A dictionary containing the optimal values for x, y, λ that maximize/minimize `f` subject to `g`.

# Examples
```julia
using PythonCall
sp = pyimport("sympy")

x, y = sp.symbols("x y")
f = x^2 + y^2
g = x + 2*y - 4
lagrangian(f, g, x, y)
```
"""
function lagrangian(f, g, x, y)
	λ = sp.symbols("λ")
	L = f - λ*g
	
	∂Lx = sp.diff(L, x)
	∂Ly = sp.diff(L, y)
	∂Lλ = sp.diff(L, λ)

	sys = [
		sp.Eq(∂Lx, 0)
		sp.Eq(∂Ly, 0)
		sp.Eq(∂Lλ, 0)
	]

	return sp.solve(sys, (x, y, λ))
end

export lagrangian